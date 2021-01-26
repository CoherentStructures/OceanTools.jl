# OceanTools.jl

[![build status][build-img]][build-url] [![coverage][codecov-img]][codecov-url]

This package was previously known as `CopernicusUtils.jl`.
Utilities for working with oceanographic datasets from the [Copernicus](http://marine.copernicus.eu/)
product family. This package can be used (amongst other things) to generate large test
cases for use with [CoherentStructures.jl](https://github.com/CoherentStructures/CoherentStructures.jl).

## Main features

Ability to load velocity fields of arbitrary space-time cubes from Copernicus CMES data files.

Fast, allocation free interpolation on regular grids supporting periodic boundaries.

## Documentation

[![stable docs][docs-stable-img]][docs-stable-url] [![dev docs][docs-dev-img]][docs-dev-url]

## Example Pictures

FTLE field for a 90 day period off the coast of Japan.

![FTLE field off the coast of Japan](https://raw.githubusercontent.com/CoherentStructures/OceanTools.jl/master/examples/ftle_plot.jpg)

Calculating ``material barriers'' on freely choosable domain, more details in the
corresponding [documentation page](https://coherentstructures.github.io/OceanTools.jl/dev/example/):

![material barriers](https://github.com/natschil/misc/raw/master/images/oceantools3.png)

## Misc

This package is in no way affiliated with the Copernicus Marine Environment Monitoring Service.

[build-img]: https://github.com/CoherentStructures/OceanTools.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/CoherentStructures/OceanTools.jl/actions?query=workflow%3ACI+branch%3Amaster

[codecov-img]: http://codecov.io/github/CoherentStructures/OceanTools.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/CoherentStructures/OceanTools.jl?branch=master

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: http://coherentstructures.github.io/OceanTools.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://coherentstructures.github.io/OceanTools.jl/stable
