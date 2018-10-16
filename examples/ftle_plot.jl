
using  Distributed, Tensors, StaticArrays, Statistics, Plots,BenchmarkTools
#addprocs(4)

import CoherentStructures
using CopernicusUtils

print("Reading in Copernicus data...\n")
ww_ocean_data = "/media/public/schillna/ocean1/worldwide_ocean_data/"

p,times = getP(ww_ocean_data, ndays=90, sshs=true)


print("Done reading in data")
@everywhere GC.gc()

LL = [100.0,15.0]
UR = [180.0,48.0]

u = @SVector [UR[1],UR[2]]
print(uv_trilinear(u,p,times[10]))
@btime uv_trilinear($u,$p,$times[10])

print(uv_tricubic(u,p,times[10]))
print(ssh_tricubic(u,p,times[10]))

#Currently takes about 9.8 microseconds returns [0.0347813, 0.0222791]
@btime uv_tricubic($u,$p,$times[10])

xs = range(130,stop=150, length=500)
ys = range(20,stop=30,length=500)
#zs = @time [fast_trilinear_earth_interpolate((@SVector [x,y]),p,times[10])[1] for x in xs, y in ys]

zs = @time [ssh_tricubic((@SVector [x,y]),p,times[1])
     for y in ys, x in xs]

zs = @time [sshVelocityRHS((@SVector [x,y]),p,times[1])[2]
     for y in ys, x in xs]
Plots.contour(xs,ys,zs)

zs = @time [uv_trilinear((@SVector [x,y]),p,times[end])[2]
     for y in ys, x in xs]
Plots.contour(xs,ys,zs)


zs = @time [uv_tricubic((@SVector [x,y]),p,times[end])[1,1]
     for y in ys, x in xs]


zs = @time [uv_tricubic((@SVector [x,y]),p,times[1])[2]
     for y in ys, x in xs]
Plots.contour(xs,ys,zs)


#xs = range(130,stop=150, length=100)
#ys = range(20,stop=30,length=100)
xs = range(LL[1],stop=UR[1],length=200)
ys = range(LL[2],stop=UR[2],length=200)
#f = x -> norm(CoherentStructures.linearized_flow(uv_tricubic,x,[times[1],times[end]],1e-8; p=p))
f = x -> maximum(eigvals(eigen(dott(
    CoherentStructures.linearized_flow(
    uv_trilinear,x,[times[1],times[end]],1e-8;
    p=p,tolerance=1e-8)[end]))))
zs = [log2(f(@SVector [x,y])) for y in ys, x in xs]
Plots.heatmap(xs,ys,zs,color=:rainbow,aspect_ratio=1)
Plots.pdf("/tmp/myoutput.pdf")





res = CoherentStructures.plot_ftle(
       uv_trilinear,p,
       [times[1],times[end]],
       LL,UR,5000,2500;
       tolerance=1e-6,
       aspect_ratio=1,
       color=:rainbow
       )
res

Plots.pdf(res,"/tmp/output_trilinear.pdf")
