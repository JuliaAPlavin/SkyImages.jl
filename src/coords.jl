abstract type AbstractProjectedCoords <: AbstractSkyCoords end

# allow projected coords with specified x/y axis direction (e.g. to make jet horizontal) ?

struct ProjectedCoords{TC,T} <: AbstractProjectedCoords
    center::TC
    xy::NTuple{2,T}
end

ConstructionBase.constructorof(::Type{<:ProjectedCoords{Val{C}}}) where {C} = (xy...,) -> ProjectedCoords(Val(C), xy)

parent_coords_type(::Type{<:ProjectedCoords}) = Any
parent_coords_type(::Type{<:ProjectedCoords{TC}}) where {TC <: AbstractSkyCoords} = TC
parent_coords_type(::Type{<:ProjectedCoords{Val{center}}}) where {center} = center isa Type ? center : typeof(center)
origin(c::ProjectedCoords) = c.center
origin(c::ProjectedCoords{Val{center}}) where {center} = center

const ProjectedCoordsS{C,T} = ProjectedCoords{Val{C}, T}


Base.convert(::Type{T}, c::T) where {T<:AbstractProjectedCoords} = c
# Base.convert(::Type{ProjectedCoords}, c::ProjectedCoordsS) = ProjectedCoords(origin(c), c.xy)
Base.convert(::Type{ProjectedCoords{Val}}, c::ProjectedCoords) = ProjectedCoords(Val(origin(c)), c.xy)
Base.convert(::Type{ProjectedCoordsS}, c::ProjectedCoords) = ProjectedCoords(Val(origin(c)), c.xy)

function Base.convert(::Type{TCto}, c::AbstractProjectedCoords) where {TCto<:AbstractSkyCoords}
    res = origin(c)
    res = @set lat(res) += c.xy[2]
    res = @set lon(res) += c.xy[1] / cos(lat(origin(c)))
    convert(TCto, res)
end


project(::Val{center}, c::AbstractSkyCoords) where {center} = ProjectedCoords(Val(center), project(center, c).xy)

function project(center::AbstractSkyCoords, c::AbstractSkyCoords)
    cc = convert(typeof(center), c)
    Δlon = Circular.center_angle(lon(cc) - lon(center))
    xy = (Δlon * cos(lat(center)), lat(cc) - lat(center))
    ProjectedCoords(center, xy)
end

SkyCoords.lon(c::AbstractProjectedCoords) = lon(origin(c)) + c.xy[1] / cos(lat(origin(c)))
SkyCoords.lat(c::AbstractProjectedCoords) = lat(origin(c)) + c.xy[2]

Base.isapprox(a::ProjectedCoords, b::ProjectedCoords; kwargs...) = isapprox(SkyCoords.SVector(lon(a), lat(a)), SkyCoords.SVector(lon(b), lat(b)); kwargs...)
