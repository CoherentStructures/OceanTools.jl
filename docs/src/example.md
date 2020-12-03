# Example

In this example, we use the [CoherentStructures.jl](https://github.com/CoherentStructures/CoherentStructures.jl)
package to work with some data derived from satelite altimetry. As is usual for
`Julia`, this will take much longer than normal to run the first time as things
are being compiled.

For the sake of simplicity we will restrict ourselves to the case of only a single
worker, but parallelising is straightforward.

## Getting data

Download velocity fields from [Copernicus CMES](https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046).
After having registered an account click on "Download", then "Download options"
and then use the "FTP Access" to download the files you want.

In our case, we have stored the files in `/foo/ww_ocean_data`, an example filename
is `dt_global_allsat_phy_l4_20070414_20190101.nc`.

!!! tip
    If you are using NetCDF files where the date is not in the filename the same way,
    adjust the `filename_match_to_date` parameter.

## Loading the data

Load the packages used:

```julia
using OceanTools, CoherentStructures, Dates, StaticArrays
```

The space-time window we are interested in:

```julia
start_date = DateTime("2006-11-25T00:00:00")
end_date = DateTime("2006-12-30T00:00:00")
LL_space = [-50.0, -50.0]
UR_space = [50.0, 50.0]
```

Read in velocity files (from "ugos" and "vgos" attributes in the NetCDF files)

```julia
schema = r"^dt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$"
ww_ocean_data = "/foo/ww_ocean_data/"

p, ust1, (Lon, Lat, times)  = read_ocean_velocities(ww_ocean_data,start_date, end_date; schema=schema, LL_space=LL_space, UR_space=UR_space)
```

Plot a heatmap of the longitudinal component of a tricubic interpolation of the velocity field.

```julia
using Plots
xs = range(-20, stop=15, length=400)
ys = range(-10, stop=5, length=400)
t = times[10]
Plots.heatmap(xs, ys, (x,y) -> uv_tricubic((@SVector [x,y]), p, t)[1], title="Longitudinal velocity [deg/day]", color=:viridis, aspect_ratio=1.0)
```

![](https://github.com/natschil/misc/raw/master/images/oceantools1.png)

Notice how areas where the velocity field is missing (i.e. on land), the velocity is zero.
Because knowing land areas is sometimes useful, the `ust1` variable contains a single time
slice of the velocity where missing values are `nan`.

The velocity fields used are derived from sea surface height measurements.
Load the sea surface height mesurements directly:

```julia
p2, _  = read_ocean_scalars(ww_ocean_data, start_date, end_date; schema=schema, LL_space=LL_space, UR_space=UR_space, scalar_field_name="sla")
Plots.heatmap(xs, ys, (x,y) -> scalar_tricubic((@SVector [x,y]), p2, t), title="Sea surface height anomaly [m]", color=:viridis, aspect_ratio=1.0)
```

![](https://github.com/natschil/misc/raw/master/images/oceantools2.png)

It is straightforward to plot trajectories of single particles, we'll make an
animation (download this [script](https://coherentstructures.github.io/CoherentStructures.jl/stable/videos/)
into "animations.jl" for the `animatemp4` function):

```julia
include("animations.jl")
x0 = [0.0,0.0]
trajectory = flow(uv_tricubic,x0,times; p=p)

frames = []
for (i,t) in enumerate(times)
    fig = Plots.heatmap(xs,ys,(x,y) -> scalar_tricubic((@SVector [x,y]),p2,t ),title="Sea surface height anomaly [m]",color=:viridis,aspect_ratio=1.0,clim=(-0.1,0.1))
    Plots.scatter!(fig,(trajectory[i][1],trajectory[i][2]),label="simulated drifter position")
    push!(frames,fig)
end
animatemp4(frames)
```

```@raw html
 <video controls="" height="100%" width="100%">
  <source src="https://github.com/natschil/misc/raw/master/videos/oceantools1.mp4" type="video/mp4" />
 Your browser does not support the video tag.
 </video>
```

The `flow` function used above came from the `CoherentStructures.jl` package
(which in turn calls `DifferentialEquations.jl`). `CoherentStructures.jl` is also
capable of fancier analyis:

```julia
vortices, singularities, bg = materialbarriers(
       uv_tricubic, range(-5,7.5,length=300), range(-40,stop=-28,length=300), range(times[2],stop=times[2]+30,length=30),
       LCSParameters(boxradius=4, indexradius=0.25, pmax=1.4,
                     merge_heuristics=[combine_20_aggressive]),
       p=p, on_torus=false);

fig = plot_vortices(vortices, singularities, (-5, -40), (7.5, -28);
    bg=bg, title="DBS field and transport barriers", showlabel=false)
```

![](https://github.com/natschil/misc/raw/master/images/oceantools3.png)

Computations on larger domains are of course possible. [Here](https://smai-jcm.centre-mersenne.org/item/SMAI-JCM_2020__6__101_0/)
is a paper that computes similar structures as above using `OceanTools.jl` and
`CoherentStructures.jl` on a global scale.
