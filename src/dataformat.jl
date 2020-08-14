

"""
    loadEarthMap(nat_earth_path)

Simple function that loads an image (from https://www.naturalearthdata.com) into a suitable array
with corresponding longitude and latitude coordinate vectors (each 1d).
Returns (array_of_pixel_values,lon,lat).
"""
function loadEarthMap(nat_earth_path)
    natearth = FileIO.load(nat_earth_path)
    ULpixel = [-179.98888888888889, 89.98888888888889]
    pixelSpacing = 0.02222222222222
    nLon = range(ULpixel[2], step=-pixelSpacing, length=size(natearth)[1])
    nLat = range(ULpixel[1], step=pixelSpacing, length=size(natearth)[2])
    return natearth, nLon, nLat
end

"""
    getLonLat(foldername,schema,[lon_label="longitude",lat_label="latitude"])

Here `foldername` is a string, and `schema` is a regular expression.
The regular expression should match filenames that should be parsed,
for example `r"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc\$"`.
Find the first file in the folder `foldername` that matches the regular expression,
and extract the longitude and latitude vectors.
Returns a vector of Float64.

"""
function getLonLat(foldername,schema,lon_label="longitude",lat_label="latitude")
    for filename in readdir(foldername)
        if match(filename, schema) != nothing
            d = NCD.Dataset(foldername * "/" * filename)
            lon = NCD.nomissing(d[lon_label][:], NaN)[:]
            lat = NCD.nomissing(d[lat_label][:], NaN)[:]
            close(d)
            return Float64.(lon), Float64.(lat)
        end
    end
    error("No file matching $schema found in $foldername")
end

"""
    _getTime(filename)

Convenience function that opens the NetCDF file `filename` and 
gets the first value of the "time" field. 
"""
function _getTime(filename)
    d = NCD.Dataset(filename)
    if !("time" ∈ d)
        error("No field \"time\" in $filename ")
    end
    t = d["time"][1]
    close(d)
    return t
end
"""
    _daysSince1950(t)
Convenience function
Takes a DateTime object `t` and calculates the number of days since 1950.
"""
function _daysSince1950(t)
    #TODO: deal with this better?
    tdays = round(t - Dates.DateTime(1950, 1, 1, 0, 0, 0), Dates.Hour).value/24.0
    return tdays
end

"""
    loadField(d,fieldname)
Returns `d[fieldname][:], _daysSince1950(d["time"][1])`
"""
function loadField(d, fieldname,date_to_int=_daysSince1950)
    t = d["time"][1]
    return d[fieldname][:], date_to_int(t)
end

"""
    _rescaleUV!(U, V, Lon, Lat, Us, Vs, numfound,remove_nan,lli,llj)

Convert lon-lat velocities from the (whole globe arrays) `U` and `V` from units m/second to degrees/day.
Save the result in (only part of globe) Us[:,:,numfound] (resp. Vs[:,:,numfound])
The integers `lli`,`llj` are (zero-based) indices of the grid square corresponding to the lower left corner of the domain.
"""
function _rescaleUV!(U, V, Lon, Lat,Us,Vs,numfound,remove_nan,lli,llj)
    R = 6371e3 # radius of the Earth in m
    s =3600*24/R
    n_small = size(Us,1)
    m_small = size(Vs,2)
    n_full = size(U,1)
    m_full = size(V,2)
    missing_value = remove_nan ? 0.0 : NaN
    
    for i in 1:n_small, j in 1:m_small
        fulli = mod1(lli+i,n_full)
        fullj = llj+j
            Us[i,j,numfound] = ismissing(U[fulli,fullj]) ? missing_value :  sec(deg2rad(Lat[fullj])) * rad2deg(U[fulli,fullj]) * s
            Vs[i,j,numfound] = ismissing(V[fulli,fullj]) ? missing_value : rad2deg(V[fulli,fullj]) * s
    end
end

function _cropScalar!(U, Lon, Lat,Us,numfound,remove_nan,lli,llj)
    n_small = size(Us,1)
    m_small = size(Vs,2)
    n_full = size(U,1)
    m_full = size(V,2)
    missing_value = remove_nan ? 0.0 : NaN
    
    for i in 1:n_small, j in 1:m_small
        fulli = mod1(lli+i,n_full)
        fullj = llj+j
            Us[i,j,numfound] = ismissing(U[fulli,fullj]) ? missing_value :  U[fulli,fullj]
    end
end


"""
    periodicGoLeft(x,start,per)

Returns a value `y` so that `x` is in [`y`, `y + per`].
And `y = start+k*per` for some integer `k`.
"""
function periodicGoLeft(x::Float64,start::Float64,per::Float64)
        if start >= x
            while start > x
                start -= per
            end
        else 
            while start < x
                start += per
            end
            start -= per
        end
        return start
end

function periodicGoRight(x::Float64,start::Float64,per::Float64)
        if start <= x
            while start < x
                start += per
            end
        else 
            while start > x
                start -= per
            end
            start += per
        end
        return start
end

function _calculateSpatialIndices(Lon,Lat,LL_space,UR_space,bounds_given)

    if Lat[1] > LL_space[2]
        minlat = Lat[1]
        error("$minlat = Lat[1] > LL_space[2]")
    end

    if Lat[end] > UR_space[2]
        maxlat = Lat[end]
        error("$maxlat = Lat[end] > UR_space[2]")
    end


    #Indices below are 0-based!

    if bounds_given
        #Let's extend Lon by 360-periodicity.
        #Which grid square is LL[1] in? Here is the index of the longitude coordinate to the left it is congruent to
        #(Note that result is 0-indexed!)
        llxi = getIndex(LL[1], Lon[1],Lon[1]+360,0.0,nx_full,periodic)[1]
        #And now the longitude coordinate of the "actual" grid square.
        llx = periodicGoLeft(LL[1], Lon[llxi+1],360.0)

        #Same thing with UR[1], but here we get the index one past to the right (we use getIndex2 instead of getIndex!).
        #Note that is possible to `wrap` around here, e.g. if UR[1] = 370...
        urxi = getIndex2(UR[1], Lon[1],Lon[1]+360,0.0,nx_full,periodic)[3]
        if (periodicGoRight(UR[1], Lon[urxi+1],360.0) - llx) > 360.0
            error("The values of UR[1] and LL[1] are too far apart!")
        end
        #Here's how many data points we need to load.
        nx_small = mod(urxi - llxi,nx_full)
        #And the x coordinate of the largest datapoint.
        urx = llx + nx_small*(Lon[2] - Lon[1])

        #As longitudes are not periodic, we don't bother periodizing here.
        llyi = getIndex(LL[2], Lat[1],Lat[1]+180,0.0,ny_full,flat)[1]
        llyi = Lat[llyi+1]
        uryi = getIndex(UR[2], Lat[1],Lat[1]+180,0.0,ny_full,flat)[1] + 2
        #In case rouding point errors pushed us off the edge:
        if UR[2] == (Lat[end] + Lat[2] - Lat[1])
            llyi = ny_full-1
            uryi = ny_full
        end
        @assert uryi <= ny_full
        ny_small = uryi - llyi
        ury = lly + ny_small*(Lat[2] - Lat[1])

        perx = 360.0
    else
        llxi = 0
        llx = Lon[1]
        #urxi = nx_full-1
        urx = Lon[1] + 360.0
        nx_small = nx_full

        llyi = 0
        lly = Lat[1]
        #uryi = ny_full-1
        ury = Lat[1] + 180.0
        ny_small = ny_full
        perx = 0.0
    end

    return llxi,llyi,llx,lly, urx,ury, nx_small,ny_small,perx,pery
end



"""
    read_ocean_velocities(foldername,start_date,end_date,boundary_t,schema; [LL_space=nothing,UR_space=nothing,...])

Reads velocity fields in a space-time window. Uses the whole globe if `LL_space` and `UR_space` are `nothing`.
Else the rectangle of points will include spatially include the rectangle bounded by `LL_space`  and `UR_space`.
The first timestep loaded is the first one in the folder `foldername` matching the regular expression `schema` where the 
"time" variable (converted to `Int` via kwarg `date_to_int`) is not less than `start_date`.
Spacing of "time" is assumed to be 1 day, change `date_to_int` if your data is different (note that the units
of the velocity fields in the files were assumed to be in m/s and are converted to deg/day ).
The range will include `end_time`
Supports 360-periodic longitude in the sense that `LL_space[1]` can larger than `UR_space[1]`.
However, it *cannot* extend by more than one period. If you are very close to one period use the whole globe
to avoid issues.
The resulting velocities are rescaled to be in degrees!
Other keyword arguments are: `lon_label` with default "longitude", `lat_label` with default "latitude",
`remove_nan` with default `true, `array_ctor` with default `SharedArray{Float64}`.
`date_to_int` with default `_daysSince1950` is the method for converting `DateTime` objects to `Int`.

Returns (res_full, Ust1, (Lon, Lat, times)) where `res_full` is the corresponding `ItpMetdata` object and Ust1 is a
slice of the x-component at the first timestep without NaNs removed
"""
function read_ocean_velocities(foldername,start_date::Int,end_date::Int,boundary_t=flat,
                                schema=r"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$",;
                                LL_space=nothing,UR_space=nothing,
                                remove_nan=true,
                                array_ctor=SharedArray{Float64},
                                lon_label="longitude",
                                lat_label="latitude",
                                date_to_int=_daysSince1950,
                                u_field_name="ugos",
                                v_field_name="vgos"
                                )
    #Get vector of Longitude and latitude coordinates
    Lon, Lat = getLonLat(ww_ocean_data,schema,lon_label=lon_label,lat_label=lat_label)
    nx_full = length(Lon)
    ny_full = length(Lat)

    @assert (LL_space == nothing) == (UR_space == nothing)
    bounds_given::Bool =  (LL_space != nothing)


    llxi,llyi,llx,lly,urx,ury, nx_small,ny_small,perx,pery = _calculateSpatialIndices(Lon,Lat,LL_space,UR_space,bounds_given)

    LL = bounds_given ? copy(LL_space) : [Lon[1], Lat[1]] 
    UR = bounds_given ? copy(UR_space) : [Lon[1] + 360, Lat[1] + 180]

    nt = round(Int,ceil(end_date - start_date))
    times = array_ctor(nt)

    LonS = array_ctor(nx_full)
    LatS = array_ctor(ny_full)
    Us = array_ctor(nx_small, ny_small, nt)
    Vs = array_ctor(nx_small, ny_small, nt)
    Ust1 = array_ctor(nx_small, ny_small, 1)

    numfound = 0
    last_time = 0
    for fname_part in readdir(ww_ocean_data)
        m = match(fname_part,schema)
        (m == nothing) && continue
        fname = ww_ocean_data * "/" * fname_part

        file_date = DateTime(m.matches[1],m.matches[2],m.matches[3])
        file_days_since_1950 = _daysSince1950(file_date)

        if file_days_since_1950 < start_date
            continue
        end
        
        if file_days_since_1950 > end_date
            break
        end

        if numfound != 0 && (last_time + 1) != file_days_since_1950
            error("Time steps are not uniform, error on file $fname_part")
        end

        numfound += 1
        d = NCD.Dataset(fname)
        U, t = loadField(d, u_field_name,date_to_int)
        V, _ = loadField(d, v_field_name,date_to_int)
        @assert t == file_days_since_1950
        _rescaleUV!(U, V, Lon, Lat,Us,Vs,numfound,remove_nan,llxi,llyi)
        if numfound == 1
            #Ust1 has the wrong type to pass to _cropScalar, therefore the roundabout way below (that is slightly inefficient).
            _cropScalar!(U, Lon, Lat,Us,numfound,false,llxi,llyi)
            Ust1 .= Us[:,:,1:1]
            _cropScalar!(U, Lon, Lat,Us,numfound,remove_nan,llxi,llyi)
        end
        close(d)
        times[numfound] = t
    end
    nt = numfound

    if last_time < end_time
        throw(BoundsError("Only read in $numfound velocities!! (last time read in $last_time < last requested time $end_time)"))
    end


    boundaryX = bounds_given ? semiperiodic : periodic
    boundaryY = outofbounds

    res_full = ItpMetadata{array_ctor}(nx_small, ny_small, nt, SVector{3}([llx, lly, times[1]]),
            SVector{3}([urx,ury, times[end] + times[2] - times[1]]),SVector{3}([perx, pery, 0.0]),
            (Us, Vs), boundaryX,boundaryY,boundary_t)


    return res_full, Ust1, (Lon, Lat, times)
end

"""
    read_ocean_scalars(args...; kwargs...)

Reads in a scalar field, otherwise like `read_ocean_velocities`. 
Resulting `data` field in the ItpMetdata is a 1-tuple with an array (of type given by `array_ctor`).
"""
function read_ocean_scalars(foldername,start_date::Int,end_date::Int,boundary_t=flat,
                                schema=r"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$",;
                                LL_space=nothing,UR_space=nothing,
                                remove_nan=true,
                                array_ctor=SharedArray{Float64},
                                lon_label="longitude",
                                lat_label="latitude",
                                date_to_int=_daysSince1950,
                                scalar_field_name="ssh",
                                )
    #Get vector of Longitude and latitude coordinates
    Lon, Lat = getLonLat(ww_ocean_data,schema,lon_label=lon_label,lat_label=lat_label)
    nx_full = length(Lon)
    ny_full = length(Lat)

    @assert (LL_space == nothing) == (UR_space == nothing)
    bounds_given::Bool =  (LL_space != nothing)


    llxi,llyi,llx,lly,urx,ury, nx_small,ny_small,perx,pery = _calculateSpatialIndices(Lon,Lat,LL_space,UR_space,bounds_given)

    LL = bounds_given ? copy(LL_space) : [Lon[1], Lat[1]] 
    UR = bounds_given ? copy(UR_space) : [Lon[1] + 360, Lat[1] + 180]

    nt = round(Int,ceil(end_date - start_date))
    times = array_ctor(nt)

    LonS = array_ctor(nx_full)
    LatS = array_ctor(ny_full)
    Us = array_ctor(nx_small, ny_small, nt)
    Ust1 = array_ctor(nx_small, ny_small, 1)
    numfound = 0
    last_time = 0
    for fname_part in readdir(ww_ocean_data)
        m = match(fname_part,schema)
        (m == nothing) && continue
        fname = ww_ocean_data * "/" * fname_part

        file_date = DateTime(m.matches[1],m.matches[2],m.matches[3])
        file_days_since_1950 = _daysSince1950(file_date)

        if file_days_since_1950 < start_date
            continue
        end
        
        if file_days_since_1950 > end_date
            break
        end

        if numfound != 0 && (last_time + 1) != file_days_since_1950
            error("Time steps are not uniform, error on file $fname_part")
        end

        numfound += 1
        d = NCD.Dataset(fname)
        U, t = loadField(d, scalar_field_name,date_to_int)
        @assert t == file_days_since_1950
        _cropScalar!(U, Lon, Lat,Us,numfound,remove_nan,llxi,llyi)
        if numfound == 1
            #Ust1 has the wrong type to pass to _cropScalar, therefore the roundabout way below (that is slightly inefficient).
            _cropScalar!(U, Lon, Lat,Us,numfound,false,llxi,llyi)
            Ust1 .= Us[:,:,1:1]
            _cropScalar!(U, Lon, Lat,Us,numfound,remove_nan,llxi,llyi)
        end
        close(d)
        times[numfound] = t
    end
    nt = numfound

    if last_time < end_time
        throw(BoundsError("Only read in $numfound velocities!! (last time read in $last_time < last requested time $end_time)"))
    end


    boundaryX = bounds_given ? semiperiodic : periodic
    boundaryY = outofbounds

    res_full = ItpMetadata{array_ctor}(nx_small, ny_small, nt, SVector{3}([llx, lly, times[1]]),
            SVector{3}([urx,ury, times[end] + times[2] - times[1]]),SVector{3}([perx, pery, 0.0]),
            (Us,), boundaryX,boundaryY,boundary_t)


    return res_full, Ust1, (Lon, Lat, times)
end

