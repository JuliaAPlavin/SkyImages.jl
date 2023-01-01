using SkyImages
using SkyImages: CoordsRectangle, project, Interp
using RectiGrids
using SkyCoords
using AxisKeys
using Test

Base.isapprox(a::T, b::T; kwargs...) where {T <: AbstractSkyCoords} = isapprox(SkyCoords.lon(a), SkyCoords.lon(b); kwargs...) && isapprox(SkyCoords.lat(a), SkyCoords.lat(b); kwargs...)

@testset "proj coords" begin
    c0 = ICRSCoords(0.1, -0.2)
    c1 = ICRSCoords(0.1 + 1e-5, -0.2 + 3e-5)
    cp = project(c0, c1)
    @test cp.center == c0
    @test cp.xy[1] ≈ 0.98 * 1e-5  rtol=1e-4
    @test cp.xy[2] ≈ 3e-5
    @test convert(ICRSCoords, cp) ≈ c1
    @test convert(GalCoords, cp) ≈ convert(GalCoords, c1)
    @test cp == cp
    @test cp ≈ cp
end

@testset "rectangle" begin
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
    @test rp.a == project(c0, r.a)
    @test rp.b == project(c0, r.b)
    gp = grid(rp; lengths=3)
    @test all(collect(gp) .≈ collect(project.(Ref(c0), grid(r; lengths=3))))

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

@testset "simple image" begin
    simg = SkyImages.load("./data/vlass.fits")
    @test size(simg) == (5329,)

    coo = ICRSCoords(3.2760228432272003, 0.21609540562015125)
    @test axiskeys(simg, :coords)[123] ≈ coo
    @test eltype(axiskeys(simg, :coords)) === ICRSCoords{Float64}
    @test simg[123] ≈ -0.001051501
    for coo in [coo, convert(GalCoords, coo)]
        @test simg(Near(coo)) ≈ -0.001051501
        @test simg(Near.([coo])) ≈ [-0.001051501]
        @test simg(Interp(coo; order=0)) ≈ -0.001051501
        @test simg(Interp.([coo]; order=0)) ≈ [-0.001051501]
        @test simg(Interp(coo; order=1)) ≈ -0.001051501  atol=1e-9
        @test simg(Interp.([coo]; order=1)) ≈ [-0.001051501]  atol=1e-9
    end

    coo = ICRSCoords(3.276022, 0.216095)
    for coo in [coo, convert(GalCoords, coo)]
        @test simg(Near(coo)) ≈ -0.001051501
        @test simg(Near.([coo])) ≈ [-0.001051501]
        @test simg(Interp(coo; order=0)) ≈ -0.001051501
        @test simg(Interp.([coo]; order=0)) ≈ [-0.001051501]
        @test simg(Interp(coo; order=1)) ≈ -0.0002471391390044121
        @test simg(Interp.([coo]; order=1)) ≈ [-0.0002471391390044121]
    end

    bbox = boundingbox(axiskeys(simg, :coords))
    @test bbox ≈ CoordsRectangle(ICRSCoords(3.2759086803666007, 0.2160905200843189), ICRSCoords(3.276266122858064, 0.21643963719077577))
    @test boundingbox(GalCoords, axiskeys(simg, :coords)).a ≈ GalCoords{Float64}(4.952041833684394, 1.2998966486360135)
    @test convert(ICRSCoords, boundingbox(ProjectedCoords, axiskeys(simg, :coords)).a) ≈ bbox.a
    @test convert(GalCoords, boundingbox(ProjectedCoords{GalCoords}, axiskeys(simg, :coords)).a) ≈ GalCoords{Float64}(4.952041833684394, 1.2998966486360135)
end


import CompatHelperLocal as CHL
CHL.@check()
