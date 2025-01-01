module MakieExt

using Makie
using SkyImages: KA_TYPE, convert_to_rectangular_KA

Makie.used_attributes(::KA_TYPE) = (:coordstype, :npix)
Makie.convert_arguments(ct::Type{<:Union{Heatmap,Image,Contour,Contourf}}, img::KA_TYPE; kwargs...) = Makie.convert_arguments(ct, convert_to_rectangular_KA(img; kwargs...))

end
