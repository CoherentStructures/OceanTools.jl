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

function read_ocean_velocities(howmany,ww_ocean_data)
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
    return Lon,Lat,Us,Vs,times
end

function read_ssh(howmany,ww_ocean_data)
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
    return Lon,Lat,sshs,times
end
