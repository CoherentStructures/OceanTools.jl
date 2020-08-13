function loadEarthMap(nat_earth_path)
    natearth = FileIO.load(nat_earth_path)
    ULpixel = [-179.98888888888889, 89.98888888888889]
    pixelSpacing = 0.02222222222222
    nLon = range(ULpixel[2], step=-pixelSpacing, length=size(natearth)[1])
    nLat = range(ULpixel[1], step=pixelSpacing, length=size(natearth)[2])
    return natearth, nLon, nLat
end

function getLonLat(foldername,schema)
    for filename in readdir(foldername)
        if match(filename, schema) != nothing
            d = NCD.Dataset(foldername * "/" * filename)
            lon = NCD.nomissing(d["longitude"][:], NaN)[:]
            lat = NCD.nomissing(d["latitude"][:], NaN)[:]
            close(d)
            return Float64.(lon), Float64.(lat)
        end
    end
    throw(AssertionError("No file matching $schema found in $foldername"))
end

function getTime(filename)
    d = NCD.Dataset(filename)
    t = d["time"][1]
    close(d)
    return t
end

function daysSince1950(t)
    #TODO: deal with this better
    tdays = round(t - Dates.DateTime(1950, 1, 1, 0, 0, 0), Dates.Hour).value/24.0
    return tdays
end

function loadField(d, filename, fieldname)
    t = d["time"][1]
    return d[fieldname][:], daysSince1950(t)
end

"""
    rescaleUV(U, V, Lon, Lat, Us, Vs, numfound)

Convert lon-lat velocities `U` and `V` from units degrees/second to kilometers/day.
Save the result in Us[:,:,numfound] (resp. Vs[:,:,numfound])
The numbers lli,llj are expected to refer to zero-based indexes!
"""
function rescaleUV(U, V, Lon, Lat,Us,Vs,numfound,remove_nan,lli,llj)
    R = 6371e3 # radius of the Earth in kilometers
    s =3600*24/R
    n_small = size(Us,1)
    m_small = size(Vs,2)
    n_full = size(U,1)
    m_full = size(V,2)
    missing_value = remove_nan ? 0.0 : NaN
    
    for i in 1:n_small, j in 1:m_small
        fulli = mod1(lli+i,n_full)
        fullj = mod1(llj+j,m_full)
            Us[i,j,numfound] = ismissing(U[fulli,fullj]) ? missing_value :  sec(deg2rad(Lat[fullj])) * rad2deg(U[fulli,fullj]) * s
            Vs[i,j,numfound] = ismissing(V[fulli,fullj]) ? missing_value : rad2deg(V[fulli,fullj]) * s
    end
end

function read_ocean_velocities(ww_ocean_data,schema,start_date,end_date,boundary;
                                LL_space=nothing,UR_space=nothing,
                                remove_nan=true,
                                array_ctor=SharedArray{Float64},
                                )
    Lon, Lat = getLonLat(ww_ocean_data,schema)
    nx_full = length(Lon)
    ny_full = length(Lat)

    @assert (LL_space == nothing) == (UR_space == nothing)
    bounds_given =  (LL_space != nothing)

    LL = bounds_given ? copy(LL_space) : [Lon[1], Lat[1]] 
    UR = bounds_given ? copy(UR_space) : [Lon[1] + 360, Lat[1] + 180]

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

    #Indices below are 0-based!

    if bounds_given
        llxi = getIndex(LL[1], Lon[1],Lon[1]+360,0.0,nx_full,periodic)[1]
        llx = periodicGoLeft(LL[1], Lon[llxi+1],360.0)

        urxi = mod(getIndex(UR[1], Lon[1],Lon[1]+360,0.0,nx_full,periodic)[1] + 2, nx_full)
        nx_small = mod(urxi - llxi,nx_full)
        urx = llx + nx_small*(Lon[2] - Lon[1])

        llyi = getIndex(LL[2], Lat[1],Lat[1]+180,0.0,ny_full,periodic)[1]
        lly = periodicGoLeft(LL[2], Lat[llyi+1],180.0)

        uryi = mod(getIndex(UR[2], Lat[1],Lat[1]+180,0.0,ny_full,periodic)[1] + 2, ny_full)
        ny_small = mod(uryi - llyi,ny_full)
        ury = lly + ny_small*(Lat[2] - Lat[1])
        perx = 360.0
        pery = 180.0
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
        pery = 0.0
    end

        
    nt = round(Int,end_date - start_date)
    times = array_ctor(nt)

    LonS = array_ctor(nx_full)
    LatS = array_ctor(ny_full)
    Us = array_ctor(nx_small, ny_small, nt)
    Vs = array_ctor(nx_small, ny_small, nt)
    Ust1 = array_ctor(nx_small, ny_small, 1)

    for fname_part in readdir(ww_ocean_data)
        m = match(fname_part,schema)
        (m == nothing) && continue
        fname = ww_ocean_data * "/" * fname_part

        file_date = DateTime(m.matches[1],m.matches[2],m.matches[3])
        file_days_since_1950 = daysSince1950(file_date)

        #TODO: do we want to introduce a safety buffer to the left/right?
        if file_days_since_1950 < start_date
            continue
        end

        if file_days_since_1950 > end_date
            break
        end

        numfound += 1
        d = NCD.Dataset(fname)
        U, t = loadField(d,fname, "ugos")
        V, _ = loadField(d,fname, "vgos")
        rescaleUV(U, V, Lon, Lat,Us,Vs,numfound,remove_nan,llxi,llyi)
        close(d)
        times[numfound] = t
    end

    if numfound != nt
        throw(BoundsError("Only read in $numfound velocities!! (required $howmany)"))
    end

    Ust1 .= Us[:,:,1:1]

    res_full = ItpMetadata(nx_small, ny_small, nt, SVector{3}([llx, lly, times[1]]),
            SVector{3}([urx,ury, times[end] + times[2] - times[1]]),SVector{3}([perx, pery, 0.0]),
            (Us, Vs), boundary[1], boundary[2], boundary[3])


        return res_full, Ust1, (Lon, Lat, times)
end

function read_ssh(howmany, ww_ocean_data,schema=r"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$";
                    remove_nan=true,
                    start_date=nothing,
                    nskip=0,
                    array_ctor=SharedArray{Float64})
    Lon, Lat = getLonLat(ww_ocean_data,schema)
=======
                    array_ctor=SharedArray{Float64})
    Lon, Lat = getLonLat(ww_ocean_data * "/" * readdir(ww_ocean_data)[1])
>>>>>>> master
    times = zeros(howmany)
    sshs = zeros(length(Lon), length(Lat), howmany)
    numfound = 0
    skipped = 0
    for fname_part in readdir(ww_ocean_data)
        fname = ww_ocean_data * "/" * fname_part
        if start_date !== nothing
            if getTime(fname) < start_date
                continue
            end
        end
        if skipped < nskip
            skipped += 1
            continue
        end
        if numfound > howmany
            break
        end
        numfound += 1
        fname = ww_ocean_data * "/" * fname_part
        d = NCD.Dataset(fname)
        ssh, t = loadField(d,fname, "adt")
        times[numfound] = t
        sshs[:,:,numfound] .= ssh
        close(d)
    end

    if numfound < howmany
        throw(BoundsError("Only read in $numfound sea surface heights!! (required $howmany)"))
    end

    nx_full = size(Lon)[1]
    ny_full = size(Lat)[1]
    nt = size(times)[1]

    LonS = array_ctor(nx_full)
    LatS = array_ctor(ny_full)
    times = array_ctor(nt)

    LonS .= Lon
    LatS .= Lat

    sshsS = array_ctor(nx_full, ny_full, nt)
    sshst1S = array_ctor(nx_full, ny_full, 1)
    sshsS .= sshs
    sshst1S .= sshs[:,:,1:1]

    if remove_nan
        for t in 1:nt, j in 1:ny_full, i in 1:nx_full
            if isnan(sshsS[i,j,t])
                sshsS[i,j,t] = 0.0
            end
        end
        #=
        for j in 1:ny_full
            for i in 1:nx_full
                if isnan(sshst1S[i,j])
                        sshsA[i,j,:] .= 0
                end
            end
        end
        =#
    end
    return LonS, LatS, sshsS, times, sshst1S
end

"""
    getP(foldername; [ndays=90, sshs=false, remove_nan=true, start_date=nothing, nskip=0, b=(periodic, periodic, flat), array_ctor=SharedArray{Float64})

Reads in ocean velocity/sea surface height data from the files in `foldername`. Files are
traversed in the order returned by `readdir`, and `ndays` files are read. If `start_date`
(type DateTime) is set, only read in files with `t` field greater than or equal
`start_date`. If `nskip` is not zero, skip `nskip` additional files. The parameter `b` is
used to define boundary behaviour in (x,y,t) direction, consult the `getIndex` function for
details. The function `array_ctor` is used to determine how the data should be stored (either
`SharedArray{Float64}`, or `zeros` for a normal array).
"""
function getPFull(foldername,start_date,end_date;
                schema=r"^nrt_global_allsat_phy_l4_([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])_.*.nc$",
                sshs=false,
                remove_nan=true,
                nskip=0,
                array_ctor=SharedArray{Float64},
                boundaryT=flat)
    if sshs
        @assert !sshs
        #TODO:fix things
        res_full = ItpMetadata(nx, ny, nt, SVector{3}([Lon[1], Lat[1], times[1]]),
            SVector{3}([360.0 + Lon[1], 180.0 + Lat[1], times[end] + times[2] - times[1]]),
            (Us, Vs), b[1], b[2], b[3])

        Lon, Lat, ssh_vals, times, sshsT1 = read_ssh(ndays, foldername,(periodic,periodic,boundaryT);
            remove_nan=remove_nan, start_date=start_date, nskip=nskip, array_ctor=array_ctor)
        Us, Vs = nothing, nothing
    else

        return read_ocean_velocities(foldername,schema,start_date=start_date, end_date=end_date, (periodic,periodic,boundaryT);
            remove_nan=remove_nan, array_ctor=array_ctor)
        ssh_vals  = nothing
    end

    if !sshs
        return res_full, Ust1, (Lon, Lat, times)
    else
        return res_full, sshst1, (Lon, Lat, times)
    end
end

"""
    getPPart(foldername,ndays,LL_space, UR_space, sshs=false,remove_nan=true,start_data=nothing,nskip=0,array_ctor=SharedArray{Float64})

Like `getPFull`, but only loads part of the dataset in space
(specifiable via lower left corner `LL_space` and upper right corner `UR_space`). 
"""
function getPPart(foldername,
                schema,
                start_date, end_date, LL_space,UR_space;
                ndays=90,
                sshs=false,
                remove_nan=true,
                nskip=0,
                array_ctor=SharedArray{Float64},
                boundaryT=flat)


    b=(semiperiodic,semiperiodic,boundaryT)

    if sshs
        @assert !sshs
        #TODO fix read_ssh so that the below works
        #res_full = ItpMetadata(nx, ny, nt, SVector{3}([Lon[1], Lat[1], times[1]]),
        #            SVector{3}([360.0 + Lon[1], 180.0 + Lat[1], times[end] + times[2] - times[1]]),
        #            ssh_vals, b[1], b[2], b[3])
        Lon, Lat, ssh_vals, times, sshsT1, res_full = read_ssh(ndays, foldername,b;
            remove_nan=remove_nan, start_date=start_date, nskip=nskip, array_ctor=array_ctor)
        Us, Vs = nothing, nothing
        return res_full, sshst1, (Lon, Lat, times)
    else
        Lon, Lat, Us, Vs, times, Ust1, res_full = read_ocean_velocities(foldername,schema,start_date=start_date, end_date=end_date, b;
            LL_space=LL_space,UR_space=UR_space,
            remove_nan=remove_nan, array_ctor=array_ctor)
        ssh_vals  = nothing
        return res_full, Ust1, (Lon, Lat, times)
    end
        
end


function restrictP(full_data, UR, LL, tspan, flow, ntraj=20, safety_factor=0.2, sshs=false)
    if !sshs
        #Calculate UR_big, LL_big
        xs = range(LL[1], stop=UR[1], length=ntraj)
        ys = range(LL[1], stop=UR[1], length=ntraj)
        grid = SVector{2}.(xs', ys)
        images = map(x -> flow(uv_trilinear, x, tspan; p=full_data), grid)
        firstcoord = x -> x[1]

        centre = 0.5*(LL + UR)
        toadd = maximum.(abs.(UR_image - centre, LL_image-centre)) * safety_factor

        LL_new = centre .- toadd
        UR_new = centre .+ toadd

        #Copy data to new SharedArray

        #make new interpolation metadata

    else
        throw(AssertionError("Not yet implemented"))
    end
end
