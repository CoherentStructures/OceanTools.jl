# OceanTools.jl
[![Build Status](https://travis-ci.org/CoherentStructures/OceanTools.jl.svg?branch=master)](https://travis-ci.org/CoherentStructures/OceanTools.jl)
[![codecov](https://codecov.io/gh/CoherentStructures/OceanTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/CoherentStructures/OceanTools.jl)

This package was previously known as `CopernicusUtils.jl`.
Utilities for working with certain [Copernicus](http://marine.copernicus.eu/) Datasets. This package can be used (amongst other things) to generate large test cases to use with [CoherentStructures.jl](https://github.com/CoherentStructures/CoherentStructures.jl).

# Main features

Ability to load velocity fields of arbitrary space-time cubes from Copernicus CMES data files. 

Fast, allocation free interpolation on regular grids supporting periodic boundaries.

# Documentation
[![][docs-dev-img]][docs-dev-url]

# Example Picture

FTLE field for a 90 day period off the coast of Japan.

<p align="center">
    <img src="https://raw.githubusercontent.com/CoherentStructures/OceanTools.jl/master/examples/ftle_plot.jpg"/>
</p>

Calculating ``material barriers'' on freely choosable domain, more details in the corresponding [documentation page](https://coherentstructures.github.io/OceanTools.jl/dev/example/) :
<p align="center">
<img src="https://github.com/natschil/misc/raw/master/images/oceantools3.png"/>
</p>

# Misc

This package is in no way affiliated with the Copernicus Marine Environment Monitoring Service.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: http://coherentstructures.github.io/OceanTools.jl/dev/





