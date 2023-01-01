struct Interp{T, O}
    val::T
    order::O
end

itp_spec(o::Int) =
    o == 0 ? BSpline(Constant()) :
    o == 1 ? BSpline(Linear()) :
    o == 2 ? BSpline(Quadratic(InPlace(OnCell()))) :
    error("unsupported order = $o")


struct WCSAxkeys{NS,NW,CT} <: AbstractVector{CT}
    wcs::WCS.WCSTransform
    size::NTuple{NS,Int}
end

WCSAxkeys(wcs, size) = WCSAxkeys{length(size),wcs.naxis,coordstype(wcs)}(wcs, size)
coordstype(w::WCSAxkeys{NS,NW,CT}) where {NS,NW,CT} = CT

Base.size(w::WCSAxkeys) = (prod(w.size),)
function Base.getindex(w::WCSAxkeys{NS,NW}, ix::Int) where {NS,NW}
    pix = (Tuple(CartesianIndices(w.size)[ix])..., ntuple(Returns(1), NW - NS)...)
    world = WCS.pix_to_world(w.wcs, collect(Float64, pix))
    wrad = NTuple{NS}(@view world[1:NS]) .|> deg2rad
    coordstype(w)(wrad...)
end

function AxisKeys.findindex(sel::Near{<:AbstractSkyCoords}, w::WCSAxkeys{NS,NW}) where {NS,NW}
    valc = convert(coordstype(w), sel.val)
    world = (SkyCoords.lon(valc) |> rad2deg, SkyCoords.lat(valc) |> rad2deg, ntuple(Returns(1), NW - 2)...)
    pix_full = WCS.world_to_pix(w.wcs, collect(Float64, world))
    pix = NTuple{NS}(@view pix_full[1:NS])
    pix_i = clamp.(round.(Int, pix), (:).(1, w.size))
    LinearIndices(w.size)[CartesianIndex(pix_i)]
end

# world_to_pix has a constant overhead of a few μs, this batch method calls it only once
function AxisKeys.findindex(sels::AbstractArray{<:Near{<:AbstractSkyCoords}}, w::WCSAxkeys{NS,NW}) where {NS,NW}
    worlds = map(sels) do sel
        valc = convert(coordstype(w), sel.val)
        world = (SkyCoords.lon(valc) |> rad2deg, SkyCoords.lat(valc) |> rad2deg, ntuple(Returns(1.0), NW - 2)...)
    end

    # world_to_pix: pass flat vector of world coords (as 2d matrix), get flat vector of pix coords
    pixs_vec_full = WCS.world_to_pix(w.wcs, collect(reinterpret(reshape, Float64, vec(worlds))))
    pixs_vec = reinterpret(reshape, NTuple{NS,Float64}, @view pixs_vec_full[1:NS, :])

    pixs = @set vec(sels) = pixs_vec
    map(pixs) do pix
        pix_i = clamp.(round.(Int, pix), (:).(1, w.size))
        LinearIndices(w.size)[CartesianIndex(pix_i)]
    end
end

AxisKeys.getkey(A, sels::AbstractArray{<:Interp{<:AbstractSkyCoords}}) = _getkey(A, sels, only(axiskeys(A)))

# world_to_pix has a constant overhead of a few μs, this batch method calls it only once
function _getkey(A, sels::AbstractArray{<:Interp{<:AbstractSkyCoords}}, w::WCSAxkeys{NS,NW}) where {NS,NW}
    worlds = map(sels) do sel
        valc = convert(coordstype(w), sel.val)
        world = (SkyCoords.lon(valc) |> rad2deg, SkyCoords.lat(valc) |> rad2deg, ntuple(Returns(1.0), NW - 2)...)
    end

    # world_to_pix: pass flat vector of world coords (as 2d matrix), get flat vector of pix coords
    pixs_vec_full = WCS.world_to_pix(w.wcs, collect(reinterpret(reshape, Float64, vec(worlds))))
    pixs_vec = reinterpret(reshape, NTuple{NS,Float64}, @view pixs_vec_full[1:NS, :])
    pixs = @set vec(sels) = pixs_vec

    @assert allequal(s.order for s in sels)
    itp = @p A |>
        reshape(__, w.size) |>
        interpolate(__, itp_spec(sels[1].order)) |>
        extrapolate(__, NaN)

    map(pixs) do pix
        itp(pix...)
    end
end


struct HealpixAxkeys{Order,TR,TT,CT} <: AbstractVector{CT}
    resolution::TR
    tree::TT

    function HealpixAxkeys(hm::Healpix.HealpixMap{<:Any,Order}; coordstype) where {Order}
        tmp = new{Order,typeof(hm.resolution),Nothing,Nothing}(hm.resolution)
        tree = nothing # NN.KDTree(map(coords_to_3d, tmp))
        return new{Order,typeof(hm.resolution),typeof(tree),coordstype}(hm.resolution, tree)
    end
end

coordstype(h::HealpixAxkeys{Order,TR,TT,CT}) where {Order,TR,TT,CT} = CT

_pix2ang(h::HealpixAxkeys{Healpix.RingOrder}, ix) = Healpix.pix2angRing(h.resolution, ix)
_pix2ang(h::HealpixAxkeys{Healpix.NestedOrder}, ix) = Healpix.pix2angNest(h.resolution, ix)
_ang2pix(h::HealpixAxkeys{Healpix.RingOrder}, theta, phi) = Healpix.ang2pixRing(h.resolution, theta, phi)
_ang2pix(h::HealpixAxkeys{Healpix.NestedOrder}, theta, phi) = Healpix.ang2pixNest(h.resolution, theta, phi)

Base.size(h::HealpixAxkeys) = (Healpix.nside2npix(h.resolution.nside),)
function Base.getindex(h::HealpixAxkeys, ix::Int)
    codec, ra = _pix2ang(h, ix)
    return coordstype(h)(ra, Healpix.colat2lat(codec))
end

function AxisKeys.findindex(sel::Near{<:AbstractSkyCoords}, h::HealpixAxkeys)
    valc = convert(coordstype(h), sel.val)
    _ang2pix(h, Healpix.lat2colat(SkyCoords.lat(valc)), mod2pi(SkyCoords.lon(valc)))
end

function AxisKeys.findindex(sel::Near{<:Union{NTuple{2},AbstractVector}}, h::HealpixAxkeys)
    @assert length(sel.val) == 2
    lon, lat = sel.val
    _ang2pix(h, Healpix.lat2colat(lat), mod2pi(lon))
end

# function Base.findall(sel::Base.Fix2{typeof(∈), <:Disk}, h::HealpixAxkeys{Healpix.RingOrder})
#     @assert length(center(sel.x)) == 2
#     NN.inrange(h.tree, coords_to_3d(center(sel.x)), AstroCatalogUtils.separation_to_chord(radius(sel.x)))
#     # Healpix.queryDiscRing(h.resolution, Healpix.lat2colat(dec), mod2pi(ra), radius(sel.x); fact=2)
# end


boundingbox(axk) = boundingbox(coordstype(axk), axk)

function boundingbox(::Type{CT}, h::HealpixAxkeys) where {CT}
    corners = (CT(0, -π / 2), CT(2π, π / 2))
    CoordsRectangle(corners...)
end

function boundingbox(::Type{CT}, w::WCSAxkeys{NS,NW}) where {NS,NW,CT<:AbstractSkyCoords}
    edgepoints = @p begin
        [
            CartesianIndex.(1, 1:w.size[2]);
            CartesianIndex.(w.size[1], 1:w.size[2]);
            CartesianIndex.(1:w.size[1], 1);
            CartesianIndex.(1:w.size[1], w.size[2]);
        ]
        map((Tuple(_)..., ntuple(Returns(1), NW - NS)...) .|> Float64)
        reinterpret(reshape, Float64, __)
        WCS.pix_to_world(w.wcs, collect(__))
        reinterpret(reshape, NTuple{NS,Float64}, @view __[1:NS, :])
        map(CT(deg2rad.(_)...))
    end
    lons = endpoints(Circular.sample_interval(lon.(edgepoints), 0..2π))
    lats = extrema(lat.(edgepoints))
    corners = CT.(lons, lats)
    rect = CoordsRectangle(corners...)
end

function boundingbox(::Type{ProjectedCoords}, axk::WCSAxkeys{NS,NW}) where {NS,NW}
    bbox = boundingbox(axk)
    origin = coordstype(axk)(deg2rad.(axk.wcs.crval[1:NS])...)
    project(origin, bbox)
end

function boundingbox(::Type{CT}, axk::WCSAxkeys{NS,NW}) where {NS,NW,CT<:ProjectedCoords}
    bbox = boundingbox(parent_coords_type(CT), axk)
    origin = convert(parent_coords_type(CT), coordstype(axk)(deg2rad.(axk.wcs.crval[1:NS])...))
    project(origin, bbox)
end
