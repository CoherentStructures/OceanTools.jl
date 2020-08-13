# Example Usage

```
using  Distributed, Tensors, StaticArrays, Statistics, Plots
addprocs(4) #Use 4 processors

import CoherentStructures
using OceanTools

#Directory containing files like nrt_global_allsat_phy_l4_20170108_20170114.nc
ww_ocean_data = "/media/public/schillna/ocean1/worldwide_ocean_data/"

#Read in data into p, missing values are set to zero
p, times = getP(ww_ocean_data, ndays=90, sshs=true, remove_nan=true)

#Garbage collection on all processors
@everywhere GC.gc()

#Bounding corners of the box on which to plot
LL = (147.0, 15.0)
UR = (180.0, 48.0)
TT = (times[1], times[end])

CoherentStructures.plot_ftle(uv_trilinear, p, TT, LL, UR, 500, 500;
       tolerance=1e-6, aspect_ratio=1, title="FTLE Field")
```
![](https://cdn.rawgit.com/CoherentStructures/OceanTools.jl/bfaf27a6/docs/example_ftle.png)
