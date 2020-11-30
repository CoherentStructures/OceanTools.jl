# Loading NetCDF datafiles

## Underlying assumptions

* There are NetCDF files in a folder (in what follows, argument `foldername`) that
  contain (time-dependent) scalar/velocity fields on a regular (in space and time)
  grid spanning the whole globe in longitude/latitude units (in particular, you do
  not care about the poles).
* Each file corresponds to a single-time snapshot of this field (and spacing between
  times are constant).
* The filenames are orded ascending with time.
* Each file has a field with longitude and latitude coordinates as a 1d array (uniform accross files). 
* You have a regular expression (e.g. `r"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$"`)
  that matches only the files you want (in what follows, argument `schema`).

## Loading velocity fields

```@docs
read_ocean_velocities
```

## Loading scalar fields

```@docs
read_ocean_scalars
```

In each case, the end result contains an `ItpMetadata` that contains all of the data needed for interpolation.

## Pseudocode example

```julia
p, _ = read_ocean_velocities(arguments)
uv_trilinear(u,p,t)
```
