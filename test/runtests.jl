using SkyImages
using SkyImages: CoordsRectangle, project
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
    g = grid(r; lengths=3)
    @test size(g) == (3, 3)
    @test g[2, 2] == ICRSCoords(0.15, -0.1)
    @test named_axiskeys(g) == (ra=[0.1, 0.15, 0.2], dec=[-0.2, -0.1, 0])

    c0 = ICRSCoords(0.11, -0.21)
    rp = project(c0, r)
    @test rp.a == project(c0, r.a)
    @test rp.b == project(c0, r.b)

    r = CoordsRectangle(GalCoords(0.1, -0.2), GalCoords(0.2, 0))
    g = grid(r; lengths=3)
    @test g[2, 2] == GalCoords(0.15, -0.1)
    @test named_axiskeys(g) == (l=[0.1, 0.15, 0.2], b=[-0.2, -0.1, 0])

    @test_throws MethodError CoordsRectangle(ICRSCoords(0.1, -0.2), GalCoords(0.2, 0))
end


import CompatHelperLocal as CHL
CHL.@check()
