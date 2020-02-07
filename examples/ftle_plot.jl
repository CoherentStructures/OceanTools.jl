using  Distributed, Tensors, StaticArrays, Statistics, Plots, BenchmarkTools
addprocs(3)

import CoherentStructures
using CopernicusUtils

print("Reading in Copernicus data...\n")
ww_ocean_data = "/media/public/schillna/ocean1/worldwide_ocean_data/"

p, Ust1, other = getP(ww_ocean_data, ndays=90, sshs=false)
times = [p.LL[3], p.UR[3]]
print("Done reading in data")
@everywhere GC.gc()
LL = (100.0, 15.0)
UR = (180.0, 48.0)

using Interpolations
const ITP = Interpolations
itptype = ITP.BSpline(ITP.Cubic(ITP.Free(ITP.OnGrid())))
const velo_1 = p.data[1]
const velo_2 = p.data[2]
@time const ui = ITP.interpolate(velo_1, itptype)
@time const vi = ITP.interpolate(velo_2, itptype)
xspan_1 = range(other[1][1], stop=other[1][end], length=length(other[1]))
yspan_1 = range(other[2][1], stop=other[2][end], length=length(other[2]))
tspan_1 = range(other[3][1], stop=other[3][end], length=length(other[3]))
UI = ITP.extrapolate(ITP.scale(ui, xspan_1, yspan_1, tspan_1), (Periodic(), Periodic(), Flat()))
VI = ITP.extrapolate(ITP.scale(vi, xspan_1, yspan_1,t span_1), (Periodic(), Periodic(), Flat()))
const p2 = (UI, VI)

@everywhere using StaticArrays
@everywhere function my_interp(u, p, t)
    v1 =  p[1](u[1], u[2], u[3], t)
    v2 =  p[2](u[1], u[2], u[3], t)
    return @SVector [v1,v2]
end

@time res = CoherentStructures.plot_ftle(
       #uv_tricubic,p,
       my_interp, p2,
       [times[1], times[end]],
       LL, UR, 100, 100;
       tolerance=1e-6,
       aspect_ratio=1,
       color=:rainbow,
       pass_on_errors=true
       );
res

Plots.pdf(res,"/tmp/output_trilinear.pdf")
