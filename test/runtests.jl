using SkyImages
using SkyImages: CoordsRectangle
using RectiGrids
using SkyCoords
using AxisKeys
using Test

@testset begin
    r = CoordsRectangle(ICRSCoords(0.1, -0.2), ICRSCoords(0.2, 0))
    g = grid(r; lengths=3)
    @test size(g) == (3, 3)
    @test g[2, 2] == ICRSCoords(0.15, -0.1)
    @test named_axiskeys(g) == (ra=[0.1, 0.15, 0.2], dec=[-0.2, -0.1, 0])

    r = CoordsRectangle(GalCoords(0.1, -0.2), GalCoords(0.2, 0))
    g = grid(r; lengths=3)
    @test g[2, 2] == GalCoords(0.15, -0.1)
    @test named_axiskeys(g) == (l=[0.1, 0.15, 0.2], b=[-0.2, -0.1, 0])

    @test_throws MethodError CoordsRectangle(ICRSCoords(0.1, -0.2), GalCoords(0.2, 0))
end


import CompatHelperLocal as CHL
CHL.@check()
