using Test, OceanTools
using Random, StaticArrays, BenchmarkTools

Random.seed!(1234)

@testset "type stability" begin
    x = rand()
    x0 = 0.0
    xf = 1.0
    nx = 123
    for boundary in instances(OceanTools.BoundaryBehaviour)
        xper = boundary == OceanTools.semiperiodic ? 3.0 : 0.0
        @inferred OceanTools.getIndex(x, x0, xf, xper, nx, boundary)
        @inferred OceanTools.getIndex2(x, x0, xf, xper, nx, boundary)
    end
    @inferred OceanTools.gooddivrem(x, nx)
end

@testset "zero allocations" begin
    xspan = range(0, stop=10.0, length=123)
    yspan = range(0, stop=10.0, length=123)
    tspan = range(0, stop=10.0, length=123)

    oob = OceanTools.outofbounds

    fu(x,y,t) = 3*x^2 + x + 2*y + π*t + 2*x*y + exp(1)*t^2  + x^2*t
    fv(x,y,t) = fu(y,x,t)

    U = [fu(x,y,t) for x in xspan, y in yspan, t in tspan]
    V = [fv(x,y,t) for x in xspan, y in yspan, t in tspan]

    metadata = @inferred OceanTools.ItpMetadata(xspan,yspan,tspan,
                  (U,V), oob, oob, oob)

    curpt = SVector{2}(10rand(2))
    t = 10rand()
    @benchmark uv_tricubic($curpt, $metadata, $t)
    @benchmark uv_trilinear($curpt, $metadata, $t)
    curmat = @SMatrix [curpt[1] 1 0; curpt[2] 0 1] 

    # type inference
    @inferred uv_trilinear(curpt, metadata, t)
    @inferred uv_tricubic(curpt, metadata, t)

    # zero allocations
    b = @benchmarkable uv_trilinear($curpt, $metadata, $t)
    r = run(b; samples=3)
    @test r.allocs == 0

    b = @benchmarkable uv_tricubic($curpt, $metadata, $t)
    r = run(b; samples=3)
    @test r.allocs == 0

    b = @benchmarkable scalar_tricubic($curpt, $metadata, $t)
    r = run(b; samples=3)
    @test r.allocs == 0

    b = @benchmarkable uv_tricubic_eqvari($curmat, $metadata, $t)
    r = run(b; samples=3)
    @test r.allocs == 0

    b = @benchmarkable scalar_tricubic_gradient($curpt, $metadata, $t)
    r = run(b; samples=3)
    @test r.allocs == 0



end
