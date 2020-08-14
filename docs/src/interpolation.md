# Interpolation on regular grids in 3D

Once you have constructed an ItpMetdata object `p` (e.g. using `read_ocean_velocities`, or by hand (see `src/interpolation.jl` for details on what the fields mean)), you can interpolate the underlying scalar/velocity field.

To do so, call `fun(u,p,t)` where `fun` is one of the functions below. Here `u` are the spatial coordinates and `t` is the time.

## Velocity fields 

```@docs
uv_tricubic
uv_trilinear
```

## Scalar fields
```@docs
scalar_trilinear
scalar_tricubic
```

## Gradients

```@docs
tricubic_gradient
```
This function has not been extensively tested/used.

## Velocity fields deriving from sea surface heights/Hamiltonians 
```@docs
ssh_rhs
```

This function has not been extensively tested/used.

## Solving the linearized flow map

```@docs
uv_tricubic_eqvari
```
This function has not been extensively tested/used.

## Boundary Behaviour

The boundary behaviour in each direction is specified using `p.boundaryX`,`p.boundaryY`, and `p.boundaryT`.
Conceptually speaking, the grid is extended to an infinite grid and interpolation is then performed on this grid.

The value `flat` means that grid points outside the original grid are assigned the value of the closest point
nt the original grid.

The value `outofbounds` means that values outside raise an error whenever they are accessed for interpolation.
Note that this means that some values that lie within the spatial bounds of the original grid will raise an error
because interpolation also requires *nearby* points.

The value `periodic` means that the grid is extended periodically in that direction. 

The value `semiperiodic` means that the grid is extended periodically (with a given periodic),
but that there can be "missing" values in between. For example, the case LL = [350,10,0] and UR = [10,20,10] 
with `p.boundaryY=semiperiodic` corresponds to the case where values are present for a small part of the globe.
Evaluating at an x coordinate of 0 (or any other multiple of 360) works (provided that the values of `y` and `t` are inbounds.
However, evaluating at an x-coordinate of 180 would raise an error.

# Interpolation on regular grids in 4D

So far, only quad-linear interpolation is supported.

```@docs
```


