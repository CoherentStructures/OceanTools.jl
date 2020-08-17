# Interpolation on regular grids in 3D

Once you have constructed an ItpMetdata object `p` (e.g. using `read_ocean_velocities`, or by hand (see below), you can interpolate the underlying scalar/velocity field.

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

There is support for calculating gradients of the interpolants in the spatial dimensions.

```@docs
scalar_tricubic_gradient
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

# Constructing an ItpMetadata object by hand

The easiest way to do this is to call the constructor `ItpMetadata(xspan,yspan,tspan, data, boundaryX,boundarY,boundaryT)`. Here `xspan`, `yspan` and `tspan` are ranges that correspond to the values of the coordinates at the datapoints. For example, if you have a $2\pi$ periodic grid with $10$ points in each direction, these should all be `range(0,stop=2\pi,length=11)[1:end-1]`.

The argument `data` is either (a) a 1-tuple containing a 3d array like object of the values at grid points or (b) a 2-tuple containing two 3d array like objects of values at grid points. Orderin should be in the form `array_like_object[x_index,y_index,t_index]` and with column-major layout. The case (a) is if you want to interpolate scalars, (b) is if you want to interpolate vectors (i.e. time-dependent spatial velocity fields)

# Interpolation on regular grids in 4D

So far, only quad-linear interpolation is supported.

```@docs
uv_quadlinear
```
This function has not been extensively tested/used.


