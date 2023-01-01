using TestItems
using TestItemRunner
@run_package_tests


@testitem "proj coords" begin
    using SkyImages: project, origin
    using SkyCoords

    Base.isapprox(a::T, b::T; kwargs...) where {T <: AbstractSkyCoords} = isapprox(SkyCoords.lon(a), SkyCoords.lon(b); kwargs...) && isapprox(SkyCoords.lat(a), SkyCoords.lat(b); kwargs...)

    c0 = ICRSCoords(0.1, -0.2)
    c1 = ICRSCoords(0.1 + 1e-5, -0.2 + 3e-5)
    cp = project(c0, c1)::ProjectedCoords
    @test origin(cp) == c0
    @test cp.xy[1] ≈ 0.98 * 1e-5  rtol=1e-4
    @test cp.xy[2] ≈ 3e-5
    @test convert(ICRSCoords, cp) ≈ c1
    @test convert(GalCoords, cp) ≈ convert(GalCoords, c1)
    @test cp == cp
    @test cp ≈ cp

    @test project(Val(c0), c1) isa ProjectedCoordsS  # deprecated
    cps = project(Val(c0), c1)::ProjectedCoords{<:Val}
    @test origin(cps) == c0
    @test cps.xy[1] ≈ 0.98 * 1e-5  rtol=1e-4
    @test cps.xy[2] ≈ 3e-5
    @test convert(ICRSCoords, cps) ≈ c1
    @test convert(GalCoords, cps) ≈ convert(GalCoords, c1)
    @test cps == cps
    @test cps ≈ cps
end

@testitem "rectangle" begin
    using SkyImages: CoordsRectangle, project
    using SkyCoords
    using RectiGrids

    Base.isapprox(a::T, b::T; kwargs...) where {T <: AbstractSkyCoords} = isapprox(SkyCoords.lon(a), SkyCoords.lon(b); kwargs...) && isapprox(SkyCoords.lat(a), SkyCoords.lat(b); kwargs...)

    r = CoordsRectangle(ICRSCoords(0.1, -0.2), ICRSCoords(0.2, 0))
    @test ICRSCoords(0.1, -0.2) ∈ r
    @test ICRSCoords(0.15, -0.15) ∈ r
    @test ICRSCoords(0.15, -0.25) ∉ r

    g = grid(r; lengths=3)
    @test size(g) == (3, 3)
    @test g[2, 2] == ICRSCoords(0.15, -0.1)
    @test named_axiskeys(g) == (ra=[0.1, 0.15, 0.2], dec=[-0.2, -0.1, 0])

    g = grid(r; lengths=(3, 5))
    @test size(g) == (3, 5)
    @test g[2, 2] == ICRSCoords(0.15, -0.15)
    @test named_axiskeys(g) == (ra = 0.1:0.05:0.2, dec = -0.2:0.05:0.0)

    c0 = ICRSCoords(0.11, -0.21)
    rp = project(c0, r)
    @test rp.a === project(c0, r.a)
    @test rp.b === project(c0, r.b)
    gp = grid(rp; lengths=3)
    @test all(collect(gp) .≈ collect(project.(Ref(c0), grid(r; lengths=3))))

    rps = project(Val(c0), r)
    @test rps.a === project(Val(c0), r.a)
    @test rps.b === project(Val(c0), r.b)
    gps = grid(rps; lengths=3)
    @test all(collect(gps) .≈ collect(project.(Val(c0), grid(r; lengths=3))))

    r = CoordsRectangle(GalCoords(0.1, -0.2), GalCoords(0.2, 0))
    g = grid(r; lengths=3)
    @test g[2, 2] == GalCoords(0.15, -0.1)
    @test named_axiskeys(g) == (l=[0.1, 0.15, 0.2], b=[-0.2, -0.1, 0])

    r = SkyImages.CoordsRectangle(ICRSCoords(5, 0.5), ICRSCoords(7, 1.3))
    g = grid(r; lengths=3)
    @test size(g) == (3, 3)
    @test g[2, 2] == ICRSCoords(6, 0.9)
    @test named_axiskeys(g) == (ra=[5, 6, 7], dec=[0.5, 0.9, 1.3])
    @test ICRSCoords(5, 1) ∈ r
    @test ICRSCoords(6, 1) ∈ r
    @test ICRSCoords(7, 1) ∈ r
    @test ICRSCoords(4, 1) ∉ r

    @test_throws MethodError CoordsRectangle(ICRSCoords(0.1, -0.2), GalCoords(0.2, 0))
end

@testitem "fits wcs image" begin
    using SkyCoords
    using SkyImages: CoordsRectangle, origin, Interp, native_rect_image
    using AxisKeys

    keyarr_equal(a, b) = a == b && named_axiskeys(a) == named_axiskeys(b)

    simg = SkyImages.load(joinpath(@__DIR__, "data/vlass.fits"))
    @test size(simg) == (5329,)

    coo = ICRSCoords(3.2760228432272003, 0.21609540562015125)
    @test axiskeys(simg, :coords)[123] ≈ coo
    @test eltype(axiskeys(simg, :coords)) === ICRSCoords{Float64}
    @test simg[123] ≈ -0.001051501
    for coo in [coo, convert(GalCoords, coo)]
        @test simg(Near(coo)) == simg[123]
        @test simg(Near.([coo])) == [simg[123]]
        @test simg(Near.([coo;;])) == [simg[123];;]
        @test keyarr_equal(simg(Near.(KeyedArray([coo;;], x=[:a], y=[5]))), KeyedArray([simg[123];;], x=[:a], y=[5]))
        @test simg(Interp(coo; order=0)) == simg[123]
        @test simg(Interp.([coo]; order=0)) == [simg[123]]
        @test simg(Interp.([coo;;]; order=0)) == [simg[123];;]
        @test keyarr_equal(simg(Interp.(KeyedArray([coo;;], x=[:a], y=[5]); order=0)), KeyedArray([simg[123];;], x=[:a], y=[5]))
        @test simg(Interp(coo; order=1)) ≈ simg[123]  atol=1e-9
    end

    coo = ICRSCoords(3.276022, 0.216095)
    for coo in [coo, convert(GalCoords, coo)]
        @test simg(Near(coo)) == simg[123]
        @test simg(Interp(coo; order=0)) == simg[123]
        @test simg(Interp(coo; order=1)) ≈ -0.0002471391390044121
    end

    @test_throws ErrorException simg(ICRSCoords(0.1, 0.2))
    @test_throws ErrorException simg(123)
    @test_throws MethodError simg(>(123))

    bbox = boundingbox(axiskeys(simg, :coords))
    @test bbox ≈ CoordsRectangle(ICRSCoords(3.2759086803666007, 0.2160905200843189), ICRSCoords(3.276266122858064, 0.21643963719077577))
    @test boundingbox(ICRSCoords, axiskeys(simg, :coords)) == bbox
    @test boundingbox(GalCoords, axiskeys(simg, :coords)).a ≈ GalCoords(4.952041833684394, 1.2998966486360135)
    @test origin(boundingbox(ProjectedCoords, axiskeys(simg, :coords)).a) == ICRSCoords(3.2754172865266287, 0.21816560027821555)
    @test origin(boundingbox(ProjectedCoords{GalCoords}, axiskeys(simg, :coords)).a) == GalCoords(4.948374097930593, 1.301730846325406)
    @test origin(boundingbox(ProjectedCoords{Val}, axiskeys(simg, :coords)).a) == ICRSCoords(3.2754172865266287, 0.21816560027821555)
    @test origin(boundingbox(ProjectedCoords{Val{GalCoords}}, axiskeys(simg, :coords)).a) == GalCoords(4.948374097930593, 1.301730846325406)
    @test convert(ICRSCoords, boundingbox(ProjectedCoords, axiskeys(simg, :coords)).a) ≈ bbox.a
    @test convert(GalCoords, boundingbox(ProjectedCoords{GalCoords}, axiskeys(simg, :coords)).a) ≈ GalCoords(4.952041833684394, 1.2998966486360135)
    @test convert(ICRSCoords, boundingbox(ProjectedCoords{Val}, axiskeys(simg, :coords)).a) ≈ bbox.a
    @test convert(GalCoords, boundingbox(ProjectedCoords{Val{GalCoords}}, axiskeys(simg, :coords)).a) ≈ GalCoords(4.952041833684394, 1.2998966486360135)

    # deprecated:
    @test origin(boundingbox(ProjectedCoordsS, axiskeys(simg, :coords)).a) == ICRSCoords(3.2754172865266287, 0.21816560027821555)
    @test origin(boundingbox(ProjectedCoordsS{GalCoords}, axiskeys(simg, :coords)).a) == GalCoords(4.948374097930593, 1.301730846325406)
    @test convert(ICRSCoords, boundingbox(ProjectedCoordsS, axiskeys(simg, :coords)).a) ≈ bbox.a
    @test convert(GalCoords, boundingbox(ProjectedCoordsS{GalCoords}, axiskeys(simg, :coords)).a) ≈ GalCoords(4.952041833684394, 1.2998966486360135)


    rimg = native_rect_image(simg)
    @test named_axiskeys(rimg).ra ≈ 3.276266446383939:-4.965847587640799e-6:3.2759089053576287
    @test named_axiskeys(rimg).dec ≈ 0.2160905962340298:4.8481455909465185e-6:0.21643966271657794
    @test rimg[12, 34] ≈ 0.05163522f0
    rimg = native_rect_image(ProjectedCoords, simg)
    @test named_axiskeys(rimg).ra ≈ 0.0008491598573101555:-4.965847587640799e-6:0.000491618831000018
    @test named_axiskeys(rimg).dec ≈ -0.0020750040441857576:4.8481455909465185e-6:-0.0017259375616376083
    @test rimg[12, 34] ≈ 0.05163522f0
end

@testitem "healpix image" begin
    using SkyCoords
    using SkyImages: CoordsRectangle, origin, Interp, native_rect_image
    using AxisKeys

    keyarr_equal(a, b) = a == b && named_axiskeys(a) == named_axiskeys(b)

    file = joinpath(@__DIR__, "data/Halpha_fwhm06_1024.fits")
    if !isfile(file)
        @warn "Cannot test healpix: image file not found"
        return
    end
    simg = SkyImages.load(file)
    @test size(simg) == (12582912,)

    coo = GalCoords(0.7892331153671623, 0.0130212012913129)
    @test axiskeys(simg, :coords)[123] ≈ coo
    @test eltype(axiskeys(simg, :coords)) === GalCoords{Float64}
    @test simg[123] ≈ 4.056696891784668
    for coo in [coo, convert(ICRSCoords, coo)]
        @test simg(Near(coo)) == simg[123]
        @test simg(Near.([coo])) == [simg[123]]
        @test simg(Near.([coo;;])) == [simg[123];;]
        @test keyarr_equal(simg(Near.(KeyedArray([coo;;], x=[:a], y=[5]))), KeyedArray([simg[123];;], x=[:a], y=[5]))
    end

    @test_throws ErrorException simg(ICRSCoords(0.1, 0.2))
    @test_throws ErrorException simg(123)
    @test_throws MethodError simg(>(123))

    bbox = boundingbox(axiskeys(simg, :coords))
    @test bbox ≈ CoordsRectangle(GalCoords(0.0, -1.5707963267948966), GalCoords(6.283185307179586, 1.5707963267948966))
    @test boundingbox(GalCoords, axiskeys(simg, :coords)) == bbox
    @test boundingbox(GalCoords, axiskeys(simg, :coords)).a ≈ GalCoords(0.0, -1.5707963267948966)
end


@testitem "" begin
    import CompatHelperLocal as CHL
    CHL.@check()
end
