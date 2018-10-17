var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#CopernicusUtils.jl-1",
    "page": "Home",
    "title": "CopernicusUtils.jl",
    "category": "section",
    "text": "Utilities for working with the Copernicus Ocean Datasets in Julia"
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "The CopernicusUtils.jl package provides a handful of helper functions to work with SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS data from Copernicus in Julia.  These are written to be as fast as possible \neffectively parallelizable\neasily usable from CoherentStructures.jlThis package was developed by Nathanael Schilling at TUM. "
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Run the following in the Julia REPL:]add https://github.com/CoherentStructures/CopernicusUtils.jl.git"
},

{
    "location": "index.html#Disclaimer-1",
    "page": "Home",
    "title": "Disclaimer",
    "category": "section",
    "text": "These functions have not been tested in detail and probably have bugs. The author of this package is in no way affiliated with Copernicus."
},

{
    "location": "index.html#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "The CopernicusUtils.jl package provides julia utilities for reading in velocity and sea surface heights (ssh) from Copernicus datasets. This functionality relies on the NetCDF.jl package. There are also functions for interpolating the resulting values (trilinear + tricubic). Tricubic interpolation is implemented using the algorithm of Lekien and Marsden\'s paper, along with a function to obtain the gradient (in space) of the interpolation function.This gives a way of approximating geostrophic sea-surface velocities with the well-known formulau = -A(y)partial_y h(xyt) v = A(y)partial_x h(xyt)where:u – Longitudinal component of the velocityv – Latitudinal component of the velocityx – Longitudey – Latitudeh – Sea-surface heightHere A(y) = gR^2 2 Omega sin y cos y  with g the gravitational constant, R the radius of the earth, Omega is the earth\'s mean angular velocity. This equation is implemented in the sshVelocityRHS function.For a list of functions that are implemented by this package, consult the exports.jl file or the source code. "
},

{
    "location": "example.html#",
    "page": "Example Usage",
    "title": "Example Usage",
    "category": "page",
    "text": ""
},

{
    "location": "example.html#Example-Usage-1",
    "page": "Example Usage",
    "title": "Example Usage",
    "category": "section",
    "text": "using  Distributed, Tensors, StaticArrays, Statistics, Plots\naddprocs(4) #Use 4 processors\n\nimport CoherentStructures\nusing CopernicusUtils\n\n#Directory containing files like nrt_global_allsat_phy_l4_20170108_20170114.nc\nww_ocean_data = \"/media/public/schillna/ocean1/worldwide_ocean_data/\"\n\n#Read in data into p, missing values are set to zero\np,times = getP(ww_ocean_data, ndays=90, sshs=true,remove_nan=true)\n\n#Garbage collection on all processors\n@everywhere GC.gc()\n\n#Bounding corners of the box on which to plot\nLL = [147.0,15.0]\nUR = [180.0,48.0]\n\nCoherentStructures.plot_ftle(\n       uv_trilinear,p,\n       [times[1],times[end]],\n       LL,UR,500,500;\n       tolerance=1e-6,\n       aspect_ratio=1,\n       title=\"FTLE Field\"\n       )(Image: )"
},

]}
