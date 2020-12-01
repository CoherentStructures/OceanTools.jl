# OceanTools.jl

*Utilities for working with the Copernicus Ocean Datasets in Julia*

## Introduction

The OceanTools.jl package provides a handful of helper functions to work with
`SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS` data from [Copernicus](http://marine.copernicus.eu/) in Julia.
These are written to be

* as fast as possible
* effectively parallelizable
* easily usable from [`CoherentStructures.jl`](https://github.com/CoherentStructures/CoherentStructures.jl)

This package was developed by Nathanael Schilling at TUM, and is maintained
by Daniel Karrasch (TUM).

## Installation

Run the following in the Julia REPL:

    ]add https://github.com/CoherentStructures/OceanTools.jl.git

## Disclaimer

These functions have not been tested in detail and probably have bugs.
The author of this package is in no way affiliated with Copernicus.

## Features

The `OceanTools.jl` package provides julia utilities for reading in velocity and
sea surface heights (ssh) from Copernicus datasets. This functionality relies on
the [`NCDatasets.jl`](https://github.com/Alexander-Barth/NCDatasets.jl) package.

There are also functions for interpolating the resulting values (trilinear + tricubic).

Tricubic interpolation is implemented using the algorithm of Lekien and Marsden's
[paper](http://www.cds.caltech.edu/~marsden/bib/2005/08-LeMa2005/LeMa2005.pdf),
along with a function to obtain the gradient (in space) of the interpolation function.

This gives a way of approximating geostrophic sea-surface velocities with the well-known formula

```math
\begin{aligned}
u &= -A(y)\partial_y h(x,y,t),\\
v &= A(y)\partial_x h(x,y,t),
\end{aligned}
```

where:

* ``u`` - longitudinal component of the velocity,
* ``v`` - latitudinal component of the velocity,
* ``x`` - longitude,
* ``y`` - latitude,
* ``h`` - sea-surface height.

Here, ``A(y) = g/(2 R \Omega \sin y)``  with $g$ the gravitational constant, ``R``
the radius of the earth, ``\Omega`` is the earth's mean angular velocity (in m/s).
This equation is implemented in the [`ssh_rhs`](@ref) function.

For a list of functions that are implemented by this package, consult the
`exports.jl` file or the source code.
