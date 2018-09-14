# CopernicusUtils.jl

*Utilities for working with the Copernicus Ocean Datasets in Julia*

# Introduction 

The CopernicusUtils.jl package provides a handful of helper functions to work with
`SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS` data from [Copernicus](http://marine.copernicus.eu/) in Julia. 
These are written to be 
   * as fast as possible 
   * effectively parallelizable
   * easily usable from [CoherentStructures.jl](https://github.com/CoherentStructures/CoherentStructures.jl)

This package was developed by Nathanael Schilling at TUM. 

# Disclaimer

These functions have not been tested in detail and probably have bugs. The author of this package is in no way affiliated with Copernicus.

# Features

The CopernicusUtils.jl package provides julia utilities for reading in velocity and sea surface heights (ssh) from Copernicus datasets.
This functionality relies on the [NetCDF.jl](https://github.com/JuliaGeo/NetCDF.jl) package. 

There are also functions for interpolating the resulting values (trilinear + tricubic). 

Tricubic interpolation is implemented using the algorithm of Lekien and Marsden's [paper](http://www.cds.caltech.edu/~marsden/bib/2005/08-LeMa2005/LeMa2005.pdf), along with a function to obtain the gradient (in space) of the interpolation function.

This gives a way of approximating geostrophic sea-surface velocities with the well-known formula

$u = -A(y)\partial_y h(x,y,t)\\ v = A(y)\partial_x h(x,y,t)$ 

where:

 $u$ -- Longitudinal component of the velocity

 $v$ -- Latitudinal component of the velocity

 $x$ -- Longitude

 $y$ -- Latitude

 $h$ -- Sea-surface height

Here $A(y) = g/R^2 2 \Omega \sin y \cos y$  with $g$ the gravitational constant, $R$ the radius of the earth, $\Omega$ is the earth's mean angular velocity. This equation is implemented in the `sshVelocityRHS` function.

For a list of functions that are implemented by this package, consult the exports.jl file or the source code. 
