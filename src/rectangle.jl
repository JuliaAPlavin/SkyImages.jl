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


function RectiGrids.grid(r::CoordsRectangle{T}; lengths) where {T}
    map(Base.splat(T),
        RectiGrids.grid(
            ra=is_lon_wrap(r) ? range(lon(r.a), lon(r.b) + 2π; length=first(lengths)) : range(lon(r.a), lon(r.b); length=first(lengths)),
            dec=range(lat(r.a), lat(r.b); length=last(lengths))
        )
    )
end

function RectiGrids.grid(r::CoordsRectangle{T}; lengths) where {T<:ProjectedCoords}
    @assert r.a.center == r.b.center
    map(xy -> ProjectedCoords(r.a.center, xy),
        RectiGrids.grid(
            Tuple,
            ra=range(r.a.xy[1], r.b.xy[1]; length=first(lengths)),
            dec=range(r.a.xy[2], r.b.xy[2]; length=last(lengths))
        )
    )
end
