struct ProjectedCoords{TC<:AbstractSkyCoords,T} <: AbstractSkyCoords
    center::TC
    xy::NTuple{2,T}
end

Base.convert(::Type{T}, c::T) where {T<:ProjectedCoords} = c

parent_coords_type(::Type{<:ProjectedCoords{TC}}) where {TC} = TC

function Base.convert(::Type{TCto}, c::ProjectedCoords) where {TCto<:AbstractSkyCoords}
    res = c.center
    res = @set lat(res) += c.xy[2]
    res = @set lon(res) += c.xy[1] / cos(lat(c.center))
    convert(TCto, res)
end

function project(center::AbstractSkyCoords, c::AbstractSkyCoords)
    cc = convert(typeof(center), c)
    xy = ((lon(cc) - lon(center)) * cos(lat(center)), lat(cc) - lat(center))
    ProjectedCoords(center, xy)
end
