function loadEarthMap(nat_earth_path)
    natearth= Images.load(nat_earth_path)
    ULpixel = [-179.98888888888889, 89.98888888888889]
    pixelSpacing = 0.02222222222222
    nLon = range(ULpixel[2], step=-1 * pixelSpacing, length=size(natearth)[1])
    nLat = range(ULpixel[1], step=pixelSpacing, length=size(natearth)[2])
    return natearth, nLon,nLat
end

function getLonLat(filename)
    return NetCDF.ncread(filename,"longitude"), NetCDF.ncread(filename,"latitude")
end

function loadUV(filename)
     U = [x == -2147483647 ? NaN : x*0.0001 for x in NetCDF.ncread(filename,"ugos")][:,:,1]
     V = [x == -2147483647 ? NaN : x*0.0001 for x in NetCDF.ncread(filename,"vgos")][:,:,1]
     t = NetCDF.ncread(filename,"time")
     return U,V,t[1]
end

function loadSSH(filename)
     ssh = [x == -2147483647 ? NaN : x*0.0001 for x in NetCDF.ncread(filename,"adt")][:,:,1]
     t = NetCDF.ncread(filename,"time")
     return ssh,t[1]
end

function rescaleUV!!(U,V,Lon,Lat)
	R = 6371e3 # Radius of the Earth in Kilometeres
	n = size(U)[1]
	m = size(U)[2]
	for j in 1:m
	    for i in 1:n
		U[i,j] = sec(deg2rad(Lat[j]))*rad2deg(U[i,j])/R*3600*24
		V[i,j] = rad2deg(V[i,j])/R*3600*24
	    end
	end
end

function read_ocean_velocities(howmany,ww_ocean_data,remove_nan=true)
    Lon, Lat = getLonLat(ww_ocean_data * "/" * readdir(ww_ocean_data)[1])
    times = zeros(howmany)
    Us = zeros(length(Lon),length(Lat),howmany)
    Vs = zeros(length(Lon),length(Lat),howmany)
    for (index,fname_part) in enumerate(readdir(ww_ocean_data))
    	if index > howmany
    	    break
    	end
    	fname = ww_ocean_data * "/" * fname_part
    	U,V,t = loadUV(fname)
    	rescaleUV!!(U,V,Lon,Lat)
    	times[index] = t
    	Us[:,:,index] .= U
    	Vs[:,:,index] .= V
    end

    sLon = size(Lon)[1]
    sLat = size(Lat)[1]
    stimes = size(times)[1]

    LonS = SharedArray{Float64}(sLon)
    LatS = SharedArray{Float64}(sLat)
    timesS = SharedArray{Float64}(stimes)

    LonS .= Lon
    LatS .= Lat
    timesS .= times

    UsS = SharedArray{Float64}(sLon,sLat,stimes)
    VsS = SharedArray{Float64}(sLon,sLat,stimes)
    Ust1S = SharedArray{Float64}(sLon,sLat,1)

    UsS .= Us
    VsS .= Vs
    Ust1S .= Us[:,:,1:1]



    UsS .= Us
    VsS .= Vs
    Ust1S .= Us[:,:,1:1]


    #Remove NaN values

    if remove_nan
        #=
        for j in 1:sLat
            for i in 1:sLon
            	if isnan(Ust1S[i,j])
                		UsS[i,j,:] .= 0
                		VsS[i,j,:] .= 0
            	end
            end
        end
        =#

        for t in 1:stimes
            for j in 1:sLat
            	for i in 1:sLon
            	    if isnan(UsS[i,j,t])
                		UsS[i,j,t] = 0.0
            	    end
            	    if isnan(Vs[i,j,t])
                		VsS[i,j,t] = 0.0
            	    end
            	end
            end
        end
    end
    return LonS,LatS,UsS,VsS,timesS,Ust1S
end

function read_ssh(howmany,ww_ocean_data,remove_nan=true)
    Lon, Lat = getLonLat(ww_ocean_data * "/" * readdir(ww_ocean_data)[1])
    times = zeros(howmany)
    sshs = zeros(length(Lon),length(Lat),howmany)
    for (index,fname_part) in enumerate(readdir(ww_ocean_data))
    	if index > howmany
    	    break
    	end
    	fname = ww_ocean_data * "/" * fname_part
    	ssh,t = loadSSH(fname)
    	times[index] = t
    	sshs[:,:,index] .= ssh
    end

    sLon = size(Lon)[1]
    sLat = size(Lat)[1]
    stimes = size(times)[1]

    LonS = SharedArray{Float64}(sLon)
    LatS = SharedArray{Float64}(sLat)
    timesS = SharedArray{Float64}(stimes)

    LonS .= Lon
    LatS .= Lat
    timesS .= times

    sshsS = SharedArray{Float64}(sLon,sLat,stimes)
    sshst1S = SharedArray{Float64}(sLon,sLat,1)
    sshsS .= sshs
    sshst1S .= sshs[:,:,1:1]


    if remove_nan
        #=
        for j in 1:sLat
            for i in 1:sLon
            	if isnan(sshst1S[i,j])
                        sshsA[i,j,:] .= 0
            	end
            end
        end
        =#

        for t in 1:stimes
            for j in 1:sLat
            	for i in 1:sLon
            	    if isnan(sshsS[i,j,t])
                		sshsS[i,j,t] = 0.0
            	    end
            	end
            end
        end
    end
    return LonS,LatS,sshsS,timesS,sshst1S
end

function getP(foldername; ndays=90, sshs=true,remove_nan=true)
    Lon,Lat,Us,Vs,times,Ust1 = read_ocean_velocities(ndays,foldername,remove_nan)
    if sshs
        Lon,Lat, sshs,times,sshsT1 = read_ssh(ndays,foldername,remove_nan)
    else
        sshs  = nothing
    end

    p = (
        Us,Vs,
        (Lon[1],Lon[end]),(Lat[1],Lat[end]), (times[1],times[end]),
        ( Ust1,Ust1,(Lon[1],Lon[end]), (Lat[1],Lat[end]),(times[1],times[2])),
        sshs
         )
    return p,times
end


