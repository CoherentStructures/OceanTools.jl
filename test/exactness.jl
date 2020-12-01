using Test, OceanTools, Random, StaticArrays, ForwardDiff

Random.seed!(1234)


@testset "exact_on_quadratic" begin
    xspan = range(0, stop=10.0, length=123)
    yspan = range(0, stop=10.0, length=123)
    tspan = range(0, stop=10.0, length=123)

    myfuncs = ((x,y,t) -> 3*x^2 + x + 2*y + π*t + 2*x*y + exp(1)*t^2  + x^2*t, (x,y,t) -> x*y*t)
    oob = OceanTools.outofbounds
    bmodes = [oob, OceanTools.periodic,OceanTools.semiperiodic, OceanTools.flat]
    for fu in myfuncs, b in bmodes
        fv(x,y,t) = fu(y,x,t)

        gu = x -> ForwardDiff.gradient(y->fu(y[1],y[2],y[3]), x)
        gv = x -> ForwardDiff.gradient(y->fv(y[1],y[2],y[3]), x)

        U = [fu(x,y,t) for x in xspan, y in yspan, t in tspan]
        V = [fv(x,y,t) for x in xspan, y in yspan, t in tspan]
        perx = (b == OceanTools.semiperiodic) ? 20 : 0.0
        metadata = @inferred OceanTools.ItpMetadata(xspan, yspan, tspan,
             (U, V), b, oob, oob; periods=(@SVector [perx, 0.0,0.0]))

        metadata2 = @inferred OceanTools.ItpMetadata(xspan, yspan, tspan,
             (U,), b, oob, oob; periods=(@SVector [perx, 0.0,0.0]))


        for i in 1:5000
            x, y, t = rand(3)

            x *= 5.0
            x += 3

            y *= 5.0
            y += 3

            t *= 5.0
            t += 3.0

            if b == OceanTools.semiperiodic
                curpt = @SVector [x + 20*(rand(Int)%5),y]
                #curpt = @SVector [x,y]
            else
                curpt = @SVector [x, y]
            end
            res2 =  uv_tricubic(curpt, metadata, t)
        
            res4 = scalar_tricubic(curpt, metadata2, t)
            @test res4 == res2[1]

            res5 = scalar_tricubic_gradient(curpt, metadata2, t)
            gradu = gu([x,y,t])
            gradv = gv([x,y,t])
            #TODO: is this sufficiently exact?

            @test res5[1] ≈ gradu[1] rtol=2e-9
            @test res5[2] ≈ gradu[2] rtol=2e-9

            curmat = @SMatrix [curpt[1] 1 0; curpt[2] 0 1] 
            res6 = uv_tricubic_eqvari(curmat, metadata, t)
            @test res6[1,1] == res2[1]
            @test res6[2,1] == res2[2]
            @test res6[1,2] ≈ res5[1] rtol=eps()
            @test res6[1,3] ≈ res5[2] rtol=eps()
            @test res6[2,2] ≈ gradv[1] rtol=2e-9
            @test res6[2,3] ≈ gradv[2] rtol=2e-9

            @test res2[1] ≈ fu(x,y,t) rtol=3e-14
            @test res2[2] ≈ fv(x,y,t) rtol=3e-14
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


    oob = OceanTools.outofbounds
    bmodes = [oob, OceanTools.periodic,OceanTools.semiperiodic, OceanTools.flat]
    for b in bmodes
        perx = (b == OceanTools.semiperiodic) ? 200.0 : 0.0

        metadata = OceanTools.ItpMetadata(length(xspan), length(yspan), length(tspan),
             (@SVector [minimum(xspan), minimum(yspan), minimum(tspan)]),
             (@SVector [maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(tspan)+step(tspan)]),(@SVector [perx,0.0,0.0]),
             (U, V), b, oob, oob)


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
end

@testset "exact_on_constant" begin
    xspan = range(0, stop=100, length=123)
    yspan = range(102.5, stop=150.0, length=15)
    tspan = range(11.5, stop=20.7, length=100)

    u, v = 5.0, -5.0
    U = fill(u, length(xspan), length(yspan), length(tspan))
    V = fill(v, length(xspan), length(yspan), length(tspan))

    for per in (OceanTools.periodic, OceanTools.flat)
        metadata = OceanTools.ItpMetadata(xspan,yspan,tspan,
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
