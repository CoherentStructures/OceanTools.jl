
using  Distributed, Tensors, StaticArrays, Statistics, Plots
addprocs(4)

import CoherentStructures
using CopernicusUtils

print("Reading in Copernicus data...\n")
ww_ocean_data = "/media/public/schillna/ocean1/worldwide_ocean_data/"

p,times = getP(ww_ocean_data, ndays=90, sshs=true)


print("Done reading in data")
@everywhere GC.gc()

LL = [100.0,15.0]
UR = [180.0,48.0]

res = CoherentStructures.plot_ftle(
       fast_trilinear_earth_interpolate,p,
       [times[1],times[end]],
       LL,UR,2000,2000;
       tolerance=1e-6,
       aspect_ratio=1,
       color=:rainbow
       )

Plots.pdf(res,"/tmp/output_trilinear.pdf")

