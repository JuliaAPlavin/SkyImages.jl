Base.@kwdef struct Interp{T, O}
    val::T
    order::O
    avoid_lons::ClosedInterval{Float64} = Inf..(-Inf)
end

Interp(val, order; kwargs...) = Interp(; val, order, kwargs...)
Interp(val; kwargs...) = Interp(; val, kwargs...)

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
coordstype(axk::WCSAxkeys{NS,NW,CT}) where {NS,NW,CT} = CT

Base.size(axk::WCSAxkeys) = (prod(axk.size),)
Base.getindex(axk::WCSAxkeys{NS,NW}, ix::Int) where {NS,NW} = @p begin
    CartesianIndices(axk.size)[ix]
    _pix_to_world(axk.wcs, __, Val(NS), Val(NW))
    coordstype(axk)(__...)
end


AxisKeys.findindex(sel, axk::WCSAxkeys) = error("Selector $sel not implemented")
AxisKeys.findindex(sel::Near{<:AbstractSkyCoords}, axk::WCSAxkeys) = AxisKeys.findindex([sel], axk) |> only

# world_to_pix has a constant overhead of a few μs, this batch method calls it only once
function AxisKeys.findindex(sels::AbstractArray{<:Near{<:AbstractSkyCoords}}, axk::WCSAxkeys{NS,NW}) where {NS,NW}
    worlds = map(sels) do sel
        valc = convert(coordstype(axk), sel.val)
        world = (SkyCoords.lon(valc) |> rad2deg, SkyCoords.lat(valc) |> rad2deg, ntuple(Returns(1.0), NW - 2)...)
    end

    # world_to_pix: pass flat vector of world coords (as 2d matrix), get flat vector of pix coords
    pixs_vec_full = WCS.world_to_pix(axk.wcs, collect(reinterpret(reshape, Float64, vec(worlds))))
    pixs_vec = reinterpret(reshape, NTuple{NS,Float64}, @view pixs_vec_full[1:NS, :])

    pixs = @set vec(sels) = pixs_vec
    map(pixs) do pix
        pix_i = clamp.(round.(Int, pix), (:).(1, axk.size))
        LinearIndices(axk.size)[CartesianIndex(pix_i)]
    end
end

_getkey(A, sels::Interp{<:AbstractSkyCoords}, axk::WCSAxkeys) = _getkey(A, [sels], axk) |> only

# world_to_pix has a constant overhead of a few μs, this batch method calls it only once
function _getkey(A, sels::AbstractArray{<:Interp{<:AbstractSkyCoords}}, axk::WCSAxkeys{NS,NW}) where {NS,NW}
    worlds = map(sels) do sel
        valc = convert(coordstype(axk), sel.val)
        valc = @modify(SkyCoords.lon(valc) |> If(∈(sel.avoid_lons))) do l
            argmin(e -> abs(e - l), endpoints(sel.avoid_lons))
        end
        world = (SkyCoords.lon(valc) |> rad2deg, SkyCoords.lat(valc) |> rad2deg, ntuple(Returns(1.0), NW - 2)...)
    end

    # world_to_pix: pass flat vector of world coords (as 2d matrix), get flat vector of pix coords
    pixs_vec_full = WCS.world_to_pix(axk.wcs, collect(reinterpret(reshape, Float64, vec(worlds))))
    pixs_vec = reinterpret(reshape, NTuple{NS,Float64}, @view pixs_vec_full[1:NS, :])
    pixs = @set vec(sels) = pixs_vec

    @assert allequal(s.order for s in sels)
    itp = @p A |>
        reshape(__, axk.size) |>
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

AxisKeys.findindex(sel, axk::HealpixAxkeys) = error("Selector $sel not implemented")
AxisKeys.findindex(sel::AbstractArray{<:Near{<:AbstractSkyCoords}}, axk::HealpixAxkeys) = AxisKeys.findindex.(sel, Ref(axk))

function AxisKeys.findindex(sel::Near{<:AbstractSkyCoords}, h::HealpixAxkeys)
    valc = convert(coordstype(h), sel.val)
    _ang2pix(h, Healpix.lat2colat(SkyCoords.lat(valc)), mod2pi(SkyCoords.lon(valc)))
end

# function Base.findall(sel::Base.Fix2{typeof(∈), <:Disk}, h::HealpixAxkeys{Healpix.RingOrder})
#     @assert length(center(sel.x)) == 2
#     NN.inrange(h.tree, coords_to_3d(center(sel.x)), AstroCatalogUtils.separation_to_chord(radius(sel.x)))
#     # Healpix.queryDiscRing(h.resolution, Healpix.lat2colat(dec), mod2pi(ra), radius(sel.x); fact=2)
# end



AxisKeys.getkey(A, sel::Interp{<:AbstractSkyCoords}) = _getkey(A, sel, only(axiskeys(A)))
AxisKeys.getkey(A, sel::AbstractArray{<:Interp{<:AbstractSkyCoords}}) = _getkey(A, sel, only(axiskeys(A)))



boundingbox(axk) = boundingbox(coordstype(axk), axk)

function boundingbox(::Type{CT}, ::HealpixAxkeys) where {CT}
    corners = (CT(0, -π / 2), CT(2π, π / 2))
    CoordsRectangle(corners...)
end

function boundingbox(::Type{CT}, axk::WCSAxkeys{NS,NW}) where {NS,NW,CT<:AbstractSkyCoords}
    edgepoints = @p begin
        [
            CartesianIndex.(1, 1:axk.size[2]);
            CartesianIndex.(axk.size[1], 1:axk.size[2]);
            CartesianIndex.(1:axk.size[1], 1);
            CartesianIndex.(1:axk.size[1], axk.size[2]);
        ]
        _pix_to_world(axk.wcs, __, Val(NS), Val(NW))
        map(coordstype(axk)(_...))
        map(convert(CT, _))
    end
    lons = endpoints(Circular.sample_interval(lon.(edgepoints), 0..2π))
    lats = extrema(lat.(edgepoints))
    corners = CT.(lons, lats)
    rect = CoordsRectangle(corners...)
end


boundingbox(::Type{ProjectedCoords}, axk::WCSAxkeys) = boundingbox(ProjectedCoords{coordstype(axk)}, axk)
boundingbox(::Type{ProjectedCoordsS}, axk::WCSAxkeys) = boundingbox(ProjectedCoordsS{coordstype(axk)}, axk)

boundingbox(::Type{CT}, axk::WCSAxkeys) where {CT <: ProjectedCoordsS} = @p begin
    boundingbox(ProjectedCoords{parent_coords_type(CT)}, axk)
    @modify(__ |> Properties()) do c
        convert(ProjectedCoordsS, c)
    end
end

function boundingbox(::Type{CT}, axk::WCSAxkeys{NS,NW}) where {NS,NW,CT<:ProjectedCoords}
    bbox = boundingbox(parent_coords_type(CT), axk)
    origin = convert(parent_coords_type(CT), coordstype(axk)(deg2rad.(axk.wcs.crval[1:NS])...))
    project(origin, bbox)
end



native_rect_image(A::KeyedArray) = _native_rect_image(coordstype(axiskeys(A, :coords).wcs), A, axiskeys(A, :coords))
native_rect_image(T::Type, A::KeyedArray) = _native_rect_image(T, A, axiskeys(A, :coords))

function _native_rect_image(::Type{T}, A, axk::WCSAxkeys{NS,NW}) where {T,NS,NW}
    let base_coord_T = T <: AbstractProjectedCoords ? parent_coords_type(T) : T
        coordstype(axk.wcs) <: base_coord_T || error("coords type $T not supported for WCS coords type $(coordstype(axk.wcs))")
    end

    data = reshape(A, axk.size)

    axworld = (
        @p( tuple.(1:axk.size[1], axk.wcs.crpix[2]) |> _pix_to_world(axk.wcs, __, Val(NS), Val(NW)) ),
        @p( tuple.(axk.wcs.crpix[1], 1:axk.size[2]) |> _pix_to_world(axk.wcs, __, Val(NS), Val(NW)) ),
    )
    axkeys = @p (
        getindex.(axworld[1], 1) |> Circular.unwrap,
        getindex.(axworld[2], 2),
    ) |> map(maybe_to_range)
    if T <: AbstractProjectedCoords
        origin = deg2rad.(axk.wcs.crval[1:NS])
        axkeys = map(.-, axkeys, origin)
    end

    names = if is_ax_separable(axk.wcs, 1) && is_ax_separable(axk.wcs, 2)
        fieldnames(coordstype(axk.wcs))
    else
        (:x, :y)
    end

    KeyedArray(data; NamedTuple{names}(axkeys)...)
end

function maybe_to_range(x::AbstractVector{<:Number}; rtol=1e-3)
    Δs = diff(x)
    abs(1 - maximum(Δs) / minimum(Δs)) < rtol ?
        range(start=first(x), stop=last(x), length=length(x)) :
        x
end


# some ad-hoc helpers
_pix_to_world(wcs::WCS.WCSTransform, pix::Union{CartesianIndex, Tuple}, args...) = _pix_to_world(wcs, [pix], args...) |> only

_pix_to_world(wcs::WCS.WCSTransform, pixs::AbstractArray{<:CartesianIndex}, args...) = @p begin
    pixs
    map(Tuple)
    _pix_to_world(wcs, __, args...)
end

_pix_to_world(wcs::WCS.WCSTransform, pixs::AbstractArray{<:Tuple{Vararg{<:Number}}}, ::Val{NS}, ::Val{NW}) where {NS,NW} = @p begin
    pixs
    map((_..., ntuple(Returns(1), NW - NS)...) .|> Float64)
    reinterpret(reshape, Float64, __)
    WCS.pix_to_world(wcs, collect(__))
    reinterpret(reshape, NTuple{NS,Float64}, @view __[1:NS, :])
    map(deg2rad.(_))
end
