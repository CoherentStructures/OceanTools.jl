using Test, CopernicusUtils, Random, StaticArrays

Random.seed!(1234)

@testset "exact_on_quadratic" begin
    xspan = range(0, stop=10.0, length=123)
    yspan = range(0, stop=10.0, length=123)
    tspan = range(0, stop=10.0, length=123)

    myfuncs = ((x,y,t) -> 3*x^2 + x + 2*y + π*t + 2*x*y + exp(1)*t^2  + x^2*t, (x,y,t) -> x*y*t)

    for fu in myfuncs
        fv(x,y,t) = fu(y,x,t)

        U = [ fu(x,y,t) for x in xspan, y in yspan, t in tspan]
        V = [ fv(x,y,t) for x in xspan, y in yspan, t in tspan]

        U = [fu(x,y,t) for x in xspan, y in yspan, t in tspan]
        V = [fv(x,y,t) for x in xspan, y in yspan, t in tspan]

        metadata = @inferred CopernicusUtils.ItpMetadata(xspan,yspan,tspan,
             (U, V), CopernicusUtils.outofbounds, CopernicusUtils.outofbounds, CopernicusUtils.outofbounds)

        for i in 1:5000
            x, y, t = rand(3)

            x *= 5.0
            x += 3

            y *= 5.0
            y += 3

            t *= 5.0
            t += 3.0

            curpt = @SVector [x,y]
            res2 =  uv_tricubic(curpt, metadata, t)

            @test res2[1] ≈ fu(x,y,t) rtol=2e-14
            @test res2[2] ≈ fv(x,y,t) rtol=2e-14
        end
    end
end


@testset "exact_on_linear" begin
    xspan = range(0, stop=100,length=123)
    yspan = range(102.5, stop=150.0,length=151)
    tspan = range(11.5, stop=20.7, length=100)

    fu(x,y,t) = x + 2*y + π*t
    fv(x,y,t) = -3*x - 5*y + 3*t

    U = [fu(x,y,t) for x in xspan, y in yspan, t in tspan]
    V = [fv(x,y,t) for x in xspan, y in yspan, t in tspan]

    metadata = CopernicusUtils.ItpMetadata(length(xspan), length(yspan), length(tspan),
         (@SVector [minimum(xspan), minimum(yspan), minimum(tspan)]),
         (@SVector [maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(tspan)+step(tspan)]),
         (U, V), 2, 2, 2)

    for i in 1:5000
        x, y, t = rand(3)

        x *= 80.0
        x += 4.0

        y *= 40
        y += 106.0

        t *= 5.0
        t += 12.0

        curpt = @SVector [x,y]
        res1 =  uv_trilinear(curpt, metadata,t)
        res2 =  uv_tricubic(curpt, metadata,t)
        @test res1[1] ≈ fu(x,y,t) rtol=2e-14
        @test res1[2] ≈ fv(x,y,t) rtol=2e-14
        @test res2[1] ≈ fu(x,y,t) rtol=2e-14
        @test res2[2] ≈ fv(x,y,t) rtol=2e-14
    end
end

@testset "exact_on_constant" begin
    xspan = range(0, stop=100, length=123)
    yspan = range(102.5, stop=150.0, length=15)
    tspan = range(11.5, stop=20.7, length=100)

    u, v = 5.0, -5.0
    U = fill(u, length(xspan), length(yspan), length(tspan))
    V = fill(v, length(xspan), length(yspan), length(tspan))

    for per in [0,1]
        metadata = CopernicusUtils.ItpMetadata(length(xspan), length(yspan), length(tspan),
             (@SVector [minimum(xspan),minimum(yspan),minimum(tspan)]),
             (@SVector [maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(tspan)+step(tspan)]),
             (U, V), per, per, per)

        for i in 1:5000
            x, y, t = rand(3)

            x *= 180.0
            x += 4.0

            y *= 140
            y += 106.0

            t *= 1005.0
            t += 12.0

            curpt = @SVector [x,y]
            res1 =  uv_trilinear(curpt, metadata,t)
            res2 =  uv_tricubic(curpt, metadata,t)
            @test res1[1] ≈ u rtol=2e-16
            @test res1[2] ≈ v rtol=2e-16

            @test res2[1] ≈ u rtol=2e-16
            @test res2[2] ≈ v rtol=2e-16
        end
    end
end
