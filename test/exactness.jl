using Random
Random.seed!(1234)
@testset "exact_on_quadratic" begin
    xspan = range(0,stop=10,length=123)
    yspan = range(0,stop=10.0,length=123)
    tspan = range(0, stop=10.0, length=123)

    fu(x,y,t) =  3*x^2 + x + 2*y + π*t + 2*x*y + exp(1)*t^2  + x^2*t
    fv(x,y,t) = fu(y,x,t)

    U = [ fu(x,y,t) for x in xspan, y in yspan, t in tspan]
    V = [ fv(x,y,t) for x in xspan, y in yspan, t in tspan]


    metadata = CopernicusUtils.ItpMetadata(length(xspan),length(yspan),length(tspan),
         (@SVector [minimum(xspan),minimum(yspan),minimum(tspan)]), (@SVector [maximum(xspan)+step(xspan),maximum(yspan)+step(yspan),maximum(tspan)+step(tspan)]), (U,V),
        2,2,2)

    ntotest = 5000

    randpts = rand(3,ntotest)

    maxe1 = 0.0
    maxe2 = 0.0

    for i in 1:ntotest
        x,y,t = randpts[:,i]

        x *= 5.0
        x += 3

        y *= 5.0
        y += 3

        t *= 5.0
        t += 3.0

        curpt = @SVector [x,y]
        res2 =  uv_tricubic(curpt, metadata,t)

        maxe1 = max(maxe1, abs(res2[1] - fu(x,y,t)))
        maxe2 = max(maxe1, abs(res2[2] - fv(x,y,t)))
    end
    @test maxe1 < 5e-9
    @test maxe2 < 5e-9
end


@testset "exact_on_linear" begin
    xspan = range(0,stop=100,length=123)
    yspan = range(102.5,stop=150.0,length=151)
    tspan = range(11.5, stop=20.7, length=100)

    fu(x,y,t) = x + 2*y + π*t
    fv(x,y,t) = -3*x - 5*y + 3*t

    U = [ fu(x,y,t) for x in xspan, y in yspan, t in tspan]
    V = [ fv(x,y,t) for x in xspan, y in yspan, t in tspan]


    metadata = CopernicusUtils.ItpMetadata(length(xspan),length(yspan),length(tspan),
         (@SVector [minimum(xspan),minimum(yspan),minimum(tspan)]), (@SVector [maximum(xspan)+step(xspan),maximum(yspan)+step(yspan),maximum(tspan)+step(tspan)]), (U,V),
        2,2,2)

    ntotest = 5000

    randpts = rand(3,ntotest)

    for i in 1:ntotest
        x,y,t = randpts[:,i]

        x *= 80.0
        x += 4.0

        y *= 40
        y += 106.0

        t *= 5.0
        t += 12.0

        curpt = @SVector [x,y]
        res1 =  uv_trilinear(curpt, metadata,t)
        res2 =  uv_tricubic(curpt, metadata,t)
        @test abs(res1[1] - fu(x,y,t)) < 5e-13
        @test abs(res1[2] - fv(x,y,t)) < 5e-13

        @test abs(res2[1] - fu(x,y,t)) < 2e-11
        @test abs(res2[2] - fv(x,y,t)) < 2e-11
    end
end

@testset "exact_on_constant" begin
    xspan = range(0,stop=100,length=123)
    yspan = range(102.5,stop=150.0,length=15)
    tspan = range(11.5, stop=20.7, length=100)

    fu(x,y,t) = 5.0
    fv(x,y,t) = -5.0

    U = [ fu(x,y,t) for x in xspan, y in yspan, t in tspan]
    V = [ fv(x,y,t) for x in xspan, y in yspan, t in tspan]


    for per in [0,1]
        metadata = CopernicusUtils.ItpMetadata(length(xspan),length(yspan),length(tspan),
             (@SVector [minimum(xspan),minimum(yspan),minimum(tspan)]), (@SVector [maximum(xspan)+step(xspan),maximum(yspan)+step(yspan),maximum(tspan)+step(tspan)]), (U,V),
            per,per,per)

        ntotest = 5000

        randpts = rand(3,ntotest)

        for i in 1:ntotest
            x,y,t = randpts[:,i]

            x *= 180.0
            x += 4.0

            y *= 140
            y += 106.0

            t *= 1005.0
            t += 12.0

            curpt = @SVector [x,y]
            res1 =  uv_trilinear(curpt, metadata,t)
            res2 =  uv_tricubic(curpt, metadata,t)
            @test abs(res1[1] - fu(x,y,t)) < 4e-13
            @test abs(res1[2] - fv(x,y,t)) < 4e-13

            @test abs(res2[1] - fu(x,y,t)) < 1e-11
            @test abs(res2[2] - fv(x,y,t)) < 1e-11
        end
    end
end
