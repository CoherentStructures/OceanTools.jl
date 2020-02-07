using Test, CopernicusUtils
using Random, StaticArrays, BenchmarkTools

@testset "zero allocations" begin
    xspan = range(0, stop=10,length=123)
    yspan = range(0, stop=10.0,length=123)
    tspan = range(0, stop=10.0, length=123)

    oob = CopernicusUtils.outofbounds

    fu(x,y,t) = 3*x^2 + x + 2*y + π*t + 2*x*y + exp(1)*t^2  + x^2*t
    fv(x,y,t) = fu(y,x,t)

    U = [fu(x,y,t) for x in xspan, y in yspan, t in tspan]
    V = [fv(x,y,t) for x in xspan, y in yspan, t in tspan]

    metadata = @inferred CopernicusUtils.ItpMetadata(
                  length(xspan), length(yspan), length(tspan),
                  (@SVector [minimum(xspan),minimum(yspan),minimum(tspan)]),
                  (@SVector [maximum(xspan)+step(xspan),maximum(yspan)+step(yspan),maximum(tspan)+step(tspan)]),
                  (U,V), oob, oob, oob)

    curpt = SVector{2}(10rand(2))
    t = 10rand()
    # @benchmark uv_tricubic($curpt, $metadata, $t)
    # @benchmark uv_trilinear($curpt, $metadata, $t)
    
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
end
