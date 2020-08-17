var documenterSearchIndex = {"docs":
[{"location":"interpolation/#Interpolation-on-regular-grids-in-3D","page":"Interpolation","title":"Interpolation on regular grids in 3D","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"Once you have constructed an ItpMetdata object p (e.g. using read_ocean_velocities, or by hand (see src/interpolation.jl for details on what the fields mean)), you can interpolate the underlying scalar/velocity field.","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"To do so, call fun(u,p,t) where fun is one of the functions below. Here u are the spatial coordinates and t is the time.","category":"page"},{"location":"interpolation/#Velocity-fields","page":"Interpolation","title":"Velocity fields","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"uv_tricubic\nuv_trilinear","category":"page"},{"location":"interpolation/#OceanTools.uv_tricubic","page":"Interpolation","title":"OceanTools.uv_tricubic","text":"uv_tricubic(x, p, t)\n\nComponent wise tricubic interpolation (Leikien-Marsden + finite differences for values not specified in their paper) of velocity field at u at time t. Velocity field stored in p.data[1] and p.data[2].\n\n\n\n\n\n","category":"function"},{"location":"interpolation/#OceanTools.uv_trilinear","page":"Interpolation","title":"OceanTools.uv_trilinear","text":"uv_trilinear(u, p, t)\n\nTrilinear interpolation of velocity field at u at time t. Velocity field stored in p.data[1] and p.data[2].\n\n\n\n\n\n","category":"function"},{"location":"interpolation/#Scalar-fields","page":"Interpolation","title":"Scalar fields","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"scalar_trilinear\nscalar_tricubic","category":"page"},{"location":"interpolation/#OceanTools.scalar_tricubic","page":"Interpolation","title":"OceanTools.scalar_tricubic","text":"scalar_tricubic(x, p, t)\n\nTricubic interpolation (Leikien-Marsden + finite differences for values not specified in their paper) of scalar field at u at time t. Scalar field stored in p.data[1]\n\n\n\n\n\n","category":"function"},{"location":"interpolation/#Gradients","page":"Interpolation","title":"Gradients","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"There is support for calculating gradients of the interpolants in the spatial dimensions.","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"scalar_tricubic_gradient","category":"page"},{"location":"interpolation/#OceanTools.scalar_tricubic_gradient","page":"Interpolation","title":"OceanTools.scalar_tricubic_gradient","text":"scalar_tricubic_gradient(u,p,t)\n\nCalculates the (spatial) gradient of the function used in scalar_tricubic\n\n\n\n\n\n","category":"function"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"This function has not been extensively tested/used.","category":"page"},{"location":"interpolation/#Velocity-fields-deriving-from-sea-surface-heights/Hamiltonians","page":"Interpolation","title":"Velocity fields deriving from sea surface heights/Hamiltonians","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"ssh_rhs","category":"page"},{"location":"interpolation/#OceanTools.ssh_rhs","page":"Interpolation","title":"OceanTools.ssh_rhs","text":"ssh_rhs(u,p,t)\n\nApproximating geostrophic sea-surface velocities with the well-known formula\n\nu = -A(y)partial_y h(xyt) v = A(y)partial_x h(xyt)\n\nwhere:\n\nu – longitudinal component of the velocity,\nv – latitudinal component of the velocity,\nx – longitude,\ny – latitude,\nh – sea-surface height.\n\nand math A(y) = g/(R^2 2 \\Omega \\sin y)\n\n\n\n\n\n","category":"function"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"This function has not been extensively tested/used.","category":"page"},{"location":"interpolation/#Solving-the-linearized-flow-map","page":"Interpolation","title":"Solving the linearized flow map","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"uv_tricubic_eqvari","category":"page"},{"location":"interpolation/#OceanTools.uv_tricubic_eqvari","page":"Interpolation","title":"OceanTools.uv_tricubic_eqvari","text":"uv_tricubic_eqvari(u,p,t)\n\nThe rhs for solving the linearized flow of the vector field (u,v) with CoherentStructures.jl\n\n\n\n\n\n","category":"function"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"This function has not been extensively tested/used.","category":"page"},{"location":"interpolation/#Boundary-Behaviour","page":"Interpolation","title":"Boundary Behaviour","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"The boundary behaviour in each direction is specified using p.boundaryX,p.boundaryY, and p.boundaryT. Conceptually speaking, the grid is extended to an infinite grid and interpolation is then performed on this grid.","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"The value flat means that grid points outside the original grid are assigned the value of the closest point nt the original grid.","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"The value outofbounds means that values outside raise an error whenever they are accessed for interpolation. Note that this means that some values that lie within the spatial bounds of the original grid will raise an error because interpolation also requires nearby points.","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"The value periodic means that the grid is extended periodically in that direction. ","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"The value semiperiodic means that the grid is extended periodically (with a given periodic), but that there can be \"missing\" values in between. For example, the case LL = [350,10,0] and UR = [10,20,10]  with p.boundaryY=semiperiodic corresponds to the case where values are present for a small part of the globe. Evaluating at an x coordinate of 0 (or any other multiple of 360) works (provided that the values of y and t are inbounds. However, evaluating at an x-coordinate of 180 would raise an error.","category":"page"},{"location":"interpolation/#Interpolation-on-regular-grids-in-4D","page":"Interpolation","title":"Interpolation on regular grids in 4D","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"So far, only quad-linear interpolation is supported.","category":"page"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"uv_quadlinear","category":"page"},{"location":"interpolation/#OceanTools.uv_quadlinear","page":"Interpolation","title":"OceanTools.uv_quadlinear","text":"uv_quadlinear(u, p, t)\n\nQuadlinear interpolation of velocity field at u at time t.\n\n\n\n\n\n","category":"function"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"This function has not been extensively tested/used.","category":"page"},{"location":"#OceanTools.jl","page":"Home","title":"OceanTools.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Utilities for working with the Copernicus Ocean Datasets in Julia","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The OceanTools.jl package provides a handful of helper functions to work with SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS data from Copernicus in Julia. These are written to be","category":"page"},{"location":"","page":"Home","title":"Home","text":"as fast as possible\neffectively parallelizable\neasily usable from CoherentStructures.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package was developed by Nathanael Schilling at TUM.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Run the following in the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"]add https://github.com/CoherentStructures/OceanTools.jl.git","category":"page"},{"location":"#Disclaimer","page":"Home","title":"Disclaimer","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"These functions have not been tested in detail and probably have bugs. The author of this package is in no way affiliated with Copernicus.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The OceanTools.jl package provides julia utilities for reading in velocity and sea surface heights (ssh) from Copernicus datasets. This functionality relies on the NetCDF.jl package.","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are also functions for interpolating the resulting values (trilinear + tricubic).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Tricubic interpolation is implemented using the algorithm of Lekien and Marsden's paper, along with a function to obtain the gradient (in space) of the interpolation function.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This gives a way of approximating geostrophic sea-surface velocities with the well-known formula","category":"page"},{"location":"","page":"Home","title":"Home","text":"u = -A(y)partial_y h(xyt) v = A(y)partial_x h(xyt)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where:","category":"page"},{"location":"","page":"Home","title":"Home","text":"u – longitudinal component of the velocity,\nv – latitudinal component of the velocity,\nx – longitude,\ny – latitude,\nh – sea-surface height.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here, A(y) = g(R 2 Omega sin y)  with g the gravitational constant, R the radius of the earth, Omega is the earth's mean angular velocity (in m/s). This equation is implemented in the ssh_rhs function.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a list of functions that are implemented by this package, consult the exports.jl file or the source code.","category":"page"},{"location":"loading_data/#Loading-NetCDF-datafiles","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"","category":"section"},{"location":"loading_data/#Underlying-Assumptions","page":"Loading NetCDF datafiles","title":"Underlying Assumptions","text":"","category":"section"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"There are NetCDF files in a folder (in what follows, argument foldername) that contain (time-dependent) scalar/velocity fields on a regular (in space and time) grid spanning the whole globe in longitude/latitude units (in particular, you do not care about the poles).","category":"page"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"Each file corresponds to a single-time snapshot of this field (and spacing between times are constant).","category":"page"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"The filenames are orded ascending with time.","category":"page"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"Each file has a field with longitude and latitude coordinates as a 1d array (uniform accross files). ","category":"page"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"You have a regular expression (e.g. r\"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$\") that matches only the files you want (in what follows, argument schema).","category":"page"},{"location":"loading_data/#Loading-velocity-fields","page":"Loading NetCDF datafiles","title":"Loading velocity fields","text":"","category":"section"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"read_ocean_velocities","category":"page"},{"location":"loading_data/#OceanTools.read_ocean_velocities","page":"Loading NetCDF datafiles","title":"OceanTools.read_ocean_velocities","text":"read_ocean_velocities(foldername,start_date,end_date,boundary_t, [schema,LL_space=nothing,UR_space=nothing,...])\n\nReads velocity fields in a space-time window. Uses the whole globe if LL_space and UR_space are nothing. Else the rectangle of points will include spatially include the rectangle bounded by LL_space  and UR_space.\n\nThe first timestep loaded is the first one in the folder foldername matching the regular expression schema where the  \"time\" variable (converted to Int via kwarg date_to_int) is not less than start_date.\n\nSpacing of \"time\" is assumed to be 1 day, change date_to_int if your data is different (note that the units of the velocity fields in the files were assumed to be in m/s and are converted to deg/day, you will have to manually rescale if spacing is not 1 day ). The range will include end_date\n\nSupports 360-periodic longitude in the sense that LL_space[1] can larger than UR_space[1]. However, it cannot extend by more than one period. If you are very close to one period use the whole globe to avoid issues.\n\nOther keyword arguments are:\n\nlon_label with default \"longitude\"\nlat_label with default \"latitude\"\nremove_nan with default true, sets missing values to 0 instead of NaN\narray_ctor with default SharedArray{Float64} to specify underlying storage to use.\ndate_to_int with default _daysSince1950 is the method for converting DateTime objects to Int\nfilename_match_to_date (if not equal to nothing) converts the match of the regex schema with the filename into the (integer) date\n\nof the contents.\n\nReturns (res_full, Ust1, (Lon, Lat, times)) where res_full is the corresponding ItpMetdata object and Ust1 is a slice of the x-component at the first timestep without NaNs removed\n\n\n\n\n\n","category":"function"},{"location":"loading_data/#Loading-scalar-fields","page":"Loading NetCDF datafiles","title":"Loading scalar fields","text":"","category":"section"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"read_ocean_scalars","category":"page"},{"location":"loading_data/#OceanTools.read_ocean_scalars","page":"Loading NetCDF datafiles","title":"OceanTools.read_ocean_scalars","text":"read_ocean_scalars(args...;scalar_field_name=\"ssh\", kwargs...)\n\nReads in a scalar field, otherwise like read_ocean_velocities.  Resulting data field in the ItpMetdata is a 1-tuple with an array (of type given by array_ctor).\n\n\n\n\n\n","category":"function"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"In each case, the end result contains an ItpMetadata that contains all of the data needed for interpolation.","category":"page"},{"location":"loading_data/#Pseudocode-example","page":"Loading NetCDF datafiles","title":"Pseudocode example","text":"","category":"section"},{"location":"loading_data/","page":"Loading NetCDF datafiles","title":"Loading NetCDF datafiles","text":"p,_ = read_ocean_velocities(arguments)\nuv_trilinear(u,p,t)","category":"page"},{"location":"example/#Example","page":"Example","title":"Example","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"In this example, we use the CoherentStructures.jl package to work with some data derived from satelite altimetry. As is usual for julia, this will take much longer than normal to run the first time  as things are being compiled.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"For the sake of simplicity we will restrict ourselves to the case of only a single worker, but parallelising is straightforward.","category":"page"},{"location":"example/#Getting-the-data","page":"Example","title":"Getting the data","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Download velocity fields from Copernicus CMES; after having made an account click on \"Download\", then \"Download Options\" and then use the \"FTP Access\" to download the files you want.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"In our case, we have stored the files in /foo/ww_ocean_data, an example filename is dt_global_allsat_phy_l4_20070414_20190101.nc.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"tip: Tip\nIf you are using NetCDF files where the date is not in the filename the same way, adjust the filename_match_to_date parameter.","category":"page"},{"location":"example/#Loading-the-data","page":"Example","title":"Loading the data","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Load the packages used:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"using OceanTools, CoherentStructures, Dates, StaticArrays","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"The space-time window we are interested in:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"start_date = DateTime(\"2006-11-25T00:00:00\")\nend_date = DateTime(\"2006-12-30T00:00:00\")\nLL_space = [-50.0,-50.0]\nUR_space = [50.0,50.0]","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Read in velocity files (from \"ugos\" and \"vgos\" attributes in the NetCDF files)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"schema=r\"^dt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$\"\nww_ocean_data=\"/foo/ww_ocean_data/\"\n\np,ust1,(Lon,Lat,times)  = read_ocean_velocities(ww_ocean_data,start_date, end_date; schema=schema,LL_space=LL_space,UR_space=UR_space)\n","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Plot a heatmap of the longitudinal component of a tricubic interpolation of the velocity field. ","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"using Plots\nxs = range(-20,stop=15,length=400)\nys = range(-10,stop=5,length=400)\nt = times[10]\nPlots.heatmap(xs,ys,(x,y) -> uv_tricubic((@SVector [x,y]),p,t )[1],title=\"Longitudinal velocity [deg/day]\",color=:viridis,aspect_ratio=1.0)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Notice how areas where the velocity field is missing (i.e. on land), the velocity is zero. Because knowing land areas is sometimes useful, the ust1 variable contains a single time slice of the velocity where missing values are nan.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"The velocity fields used are derived from sea surface height measurements. load the sea surface height mesurements directly:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"p2,_  = read_ocean_scalars(ww_ocean_data,start_date, end_date; schema=schema,LL_space=LL_space,UR_space=UR_space,scalar_field_name=\"sla\")\nPlots.heatmap(xs,ys,(x,y) -> scalar_tricubic((@SVector [x,y]),p2,t ),title=\"Sea surface height anomaly [m]\",color=:viridis,aspect_ratio=1.0)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"It is straightforward to plot trajectories of single particles, we'll make an animation (download this script into \"animations.jl\" for the animatemp4 function):","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"include(\"animations.jl\")\nx0 = [0.0,0.0]\ntrajectory = flow(uv_tricubic,x0,times; p=p)\n\nframes = []\nfor (i,t) in enumerate(times)\n    fig = Plots.heatmap(xs,ys,(x,y) -> scalar_tricubic((@SVector [x,y]),p2,t ),title=\"Sea surface height anomaly [m]\",color=:viridis,aspect_ratio=1.0,clim=(-0.1,0.1))\n    Plots.scatter!(fig,(trajectory[i][1],trajectory[i][2]),label=\"simulated drifter position\")\n    push!(frames,fig)\nend\nanimatemp4(frames)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":" <video controls=\"\" height=\"100%\" width=\"100%\">\n  <source src=\"https://github.com/natschil/misc/raw/master/videos/oceantools1.mp4\" type=\"video/mp4\" />\n Your browser does not support the video tag.\n </video>","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"The flow function used above came from the CoherentStructures.jl package (which in turn calls DifferentialEquations.jl). CoherentStructures.jl is also capable of fancier analyis:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"vortices, singularities, bg = materialbarriers(\n       uv_tricubic, range(-5,7.5,length=300), range(-40,stop=-28,length=300), range(times[2],stop=times[2]+30,length=30),\n       LCSParameters(boxradius=4, indexradius=0.25, pmax=1.4,\n                     merge_heuristics=[combine_20_aggressive]),\n       p=p, on_torus=false);\n\nfig = plot_vortices(vortices, singularities, (-5, -40), (7.5, -28);\n    bg=bg, title=\"DBS field and transport barriers\", showlabel=false)\n","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Computations on larger domains are of course possible. Here is a paper that computes similar structures as above using OceanTools.jl and CoherentStructures.jl on a global scale.","category":"page"}]
}
