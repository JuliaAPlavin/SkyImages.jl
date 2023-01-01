function coordstype(wcs::WCS.WCSTransform)
    if startswith(wcs.ctype[1], "GLON")
        @assert startswith(wcs.ctype[2], "GLAT")
        return GalCoords{Float64}
    end
    if startswith(wcs.ctype[1], "RA")
        @assert startswith(wcs.ctype[2], "DEC")
        wcs.radesys == "ICRS" && return ICRSCoords{Float64}
        wcs.radesys == "FK5" && return FK5Coords{2000,Float64}
        error("unknown radesys = $(wcs.radesys)")
    end
    error("unknown ctypes = $(wcs.ctype)")
end


is_image_hdu(hdu) = hdu isa ImageHDU && length(hdu) > 1
is_healpix_hdu(hdu) = !(hdu isa ImageHDU) && getkey(read_header(hdu), "PIXTYPE", nothing) == "HEALPIX"

function is_ax_separable(wcs::WCS.WCSTransform, i::Int)
    wa = WCS.pix_to_world(wcs, wcs.crpix)
    wb = WCS.pix_to_world(wcs, modify(x -> x + 1, wcs.crpix, @optic _[i]))
    wΔ = delete(wa .- wb, @optic _[i])
    isapprox(wΔ, zero(wΔ), atol=√eps(eltype(wΔ)))
end
