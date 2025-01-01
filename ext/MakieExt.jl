module MakieExt

using MakieExtra
using SkyImages: KA_TYPE, convert_to_rectangular_KA

Makie.used_attributes(::KA_TYPE) = (:coordstype,)
MakieExtra.@define_plotfunc_conv (heatmap, image, contour, contourf) KA_TYPE
MakieExtra._convert_arguments_singlestep(::Type{<:Union{Heatmap,Image,Contour,Contourf}}, img::KA_TYPE; kwargs...) = (convert_to_rectangular_KA(img; kwargs...),)

end
