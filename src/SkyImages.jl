module SkyImages

import ConstructionBase
using AccessorsExtra
using DataPipes
using AxisKeys
using SkyCoords
using SkyCoords: lat, lon
using FITSIO: FITS, read_header, ImageHDU
using Interpolations
using DirectionalStatistics: Circular
using IntervalSets
import WCS
import Healpix
import RectiGrids

export ProjectedCoords, ProjectedCoordsS, eval_at_coords, boundingbox


include("coords.jl")
include("rectangle.jl")
include("wcs_helpers.jl")
include("axkeys.jl")


function load(path)
    FITS(path) do f
        hdu = only(Iterators.filter(hdu -> is_image_hdu(hdu) || is_healpix_hdu(hdu), f))
        if is_healpix_hdu(hdu)
            column = 1  # XXX: hardcoded
            T = Float64  # XXX: hardcoded
            hmap = Healpix.readMapFromFITS(path, column, T)
            nu_ka = KeyedArray(hmap.pixels; coords=HealpixAxkeys(hmap; coordstype=GalCoords{Float64}))
        else
            wcs = read_header(hdu, String) |> WCS.from_header |> first
            @assert all(i -> is_ax_separable(wcs, i), 3:wcs.naxis)
            @assert all(==(1), size(hdu)[3:end])
            KeyedArray(vec(read(hdu)); coords=WCSAxkeys(wcs, size(hdu)[1:2]))
        end
    end
end

function eval_at_coords(img, coords; order=nothing, kwargs...)
    if isnothing(order)
        @assert isempty(kwargs)
        @modify(vec(coords)) do xs
            img(Near.(xs))
        end
    else
        @modify(vec(coords)) do xs
            img(Interp.(xs; order, kwargs...))
        end
    end
end



# piracy - to be upstreamed?
AxisKeys.keys_view(keys::Tuple, inds::Tuple{<:AbstractArray}) = ()
AxisKeys.keys_view(keys::Tuple, inds::Tuple{<:KeyedArray}) = axiskeys(only(inds))
# AxisKeys.keys_getindex(keys::Tuple, inds::Tuple{<:AbstractArray}) = ()
# AxisKeys.keys_getindex(keys::Tuple, inds::Tuple{<:KeyedArray}) = axiskeys(only(inds))
AxisKeys.NamedDims.remaining_dimnames_from_indexing(dimnames_::Tuple, inds::Tuple{<:KeyedArray}) = dimnames(only(inds))


end
