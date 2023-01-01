struct CoordsRectangle{T<:AbstractSkyCoords} # <: IntervalSets.Domain{T}
    a::T
    b::T
end

is_lon_wrap(r::CoordsRectangle) = lon(r.a) ≈ lon(r.b) || lon(r.a) > lon(r.b)

function Base.in(x::AbstractSkyCoords, r::CoordsRectangle{T}) where {T}
    xc = convert(T, x)
    lat(r.a) <= lat(xc) <= lat(r.b) || return false
    is_lon_wrap(r) ?
    lon(r.a) <= lon(xc) || lon(xc) <= lon(r.b) :
    lon(r.a) <= lon(xc) <= lon(r.b)
end

project(center, rect) = @modify(x -> project(center, x), rect |> Properties())

Base.isapprox(a::T, b::T; kwargs...) where {T <: CoordsRectangle} = isapprox(a.a, b.a; kwargs...) && isapprox(a.b, b.b; kwargs...)


function RectiGrids.grid(r::CoordsRectangle{T}; lengths) where {T}
    lons = range(lon(r.a), lon(r.b) + (is_lon_wrap(r) ? 2π : 0); length=first(lengths))
    lats = range(lat(r.a), lat(r.b); length=last(lengths))
    RectiGrids.grid(T; NamedTuple{fieldnames(T)}((lons, lats))...)
end

function RectiGrids.grid(r::CoordsRectangle{T}; lengths) where {T<:ProjectedCoords}
    @assert origin(r.a) === origin(r.b)
    lons = range(r.a.xy[1], r.b.xy[1]; length=first(lengths))
    lats = range(r.a.xy[2], r.b.xy[2]; length=last(lengths))
    map(xy -> ProjectedCoords(r.a.center, xy),
        RectiGrids.grid(Tuple; NamedTuple{fieldnames(parent_coords_type(T))}((lons, lats))...)
    )
end

function RectiGrids.grid(r::CoordsRectangle{T}; lengths) where {T<:ProjectedCoordsS}
    @assert origin(r.a) === origin(r.b)
    lons = range(r.a.xy[1], r.b.xy[1]; length=first(lengths))
    lats = range(r.a.xy[2], r.b.xy[2]; length=last(lengths))
    RectiGrids.grid(T; NamedTuple{fieldnames(parent_coords_type(T))}((lons, lats))...)
end
