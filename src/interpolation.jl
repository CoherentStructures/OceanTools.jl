
#Constants needed for tricubic interpolation
include("coeff.jl")

#Just divrem, but casts the first result to Int
@inline function gooddivrem(x::T,y::Int64)::Tuple{Int64,T} where T
    a,b = divrem(x,y)
    return Base.unsafe_trunc(Int,a), b
end

#=
function gooddivrem(x::ForwardDiff.Dual, y)
    a,b = divrem(x,y)
    return Int(ForwardDiff.value(a)), b
end
=#

function fast_trilinear_earth_interpolate(
        u::SArray{Tuple{2},T,1,2},p,tin::Float64
        )::SArray{Tuple{2},T,1,2} where {T<:Real}
    Us = p[1]
    Vs = p[2]
    return fast_trilinear_earth_interpolate_internal(u,p,tin,Us,Vs)
end

function fast_trilinear_earth_interpolate_internal(
        u::SArray{Tuple{2},T,1,2},p,tin::Float64,Us::U,Vs::V
        )::SArray{Tuple{2},T,1,2} where {T<:Real,U,V}
    nx::Int64 = size(Us)[1]
    ny::Int64 = size(Us)[2]
    nt::Int64 = size(Us)[3]
    #Get the spatial bounds from p
    ll1::Float64,ur1::Float64  = p[3]
    ll2::Float64,ur2::Float64 = p[4]
    t0::Float64,tf::Float64 = p[5]
    t::Float64 = tin
    if t > tf
        t = tf
    end
    if t < t0
        t = t0
    end

    xindex::Int64, xcoord::T = gooddivrem((mod((u[1] - ll1), 360)*nx)/360.0,1)
    yindex::Int64, ycoord::T = gooddivrem((mod((u[2] - ll2), 180)*ny)/180.0,1)
    #
    tindex::Int64, tcoord::Float64 = gooddivrem((nt-1)*(t-t0)/(tf-t0),1)
    tindex += 1
    #Make sure we don't go out of bounds
    tpp = tindex + 1
    if tpp > nt
        tpp = nt
    end

    #Actual interpolation for u
    #TODO: Rewrite this with the earthIndex function (see below)
    r1u::T =  Us[xindex+1,yindex + 1,tindex ]*(1 - xcoord) +
                    Us[ (xindex + 1) % nx + 1,yindex + 1, tindex ]*xcoord
    r2u::T =  Us[xindex+1,(yindex + 1) % ny + 1,tindex ]*(1 - xcoord) +
                    Us[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tindex ]*xcoord
    r3u::T =  Us[xindex+1,yindex + 1,tpp ]*(1 - xcoord) +
                    Us[ (xindex + 1) % nx + 1,yindex + 1, tpp ]*xcoord
    r4u::T =  Us[xindex+1,(yindex + 1) % ny + 1,tpp ]*(1 - xcoord) +
                    Us[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tpp ]*xcoord
    res1::T =  (
        (1-tcoord)*((1-ycoord)*r1u + ycoord*r2u)
         + tcoord*((1-ycoord)*r3u + ycoord*r4u))

    #For v
    r1v::T =  Vs[xindex+1,yindex + 1,tindex ]*(1 - xcoord) +
                    Vs[ (xindex + 1) % nx + 1,yindex + 1, tindex ]*xcoord
    r2v::T =  Vs[xindex+1,(yindex + 1) % ny + 1,tindex ]*(1 - xcoord) +
                    Vs[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tindex ]*xcoord
    r3v::T =  Vs[xindex+1,yindex + 1,tpp ]*(1 - xcoord) +
                    Vs[ (xindex + 1) % nx + 1,yindex + 1, tpp ]*xcoord
    r4v::T =  Vs[xindex+1, (yindex + 1) % ny + 1,tpp ]*(1 - xcoord) +
                    Vs[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tpp ]*xcoord
    res2::T =  (
        (1-tcoord)*((1-ycoord)*r1v + ycoord*r2v)
         + tcoord*((1-ycoord)*r3v + ycoord*r4v))
    return SArray{Tuple{2},T,1,2}((res1,res2))
end

#Periodic boundary in x and y directions, constant boundary in z direction.
@inline @inbounds function earthIndex(x::Tuple{Int64,Int64,Int64},nx::Int64,ny::Int64,nt::Int64)::CartesianIndex{3}
    return CartesianIndex( (mod(x[1],nx) + 1 ,mod(x[2],ny) + 1, max(min(x[3], nt-1 ),0) + 1))
end

#Some basic tricubic interpolation
@inbounds @inline function base_tricubic_interpolation(
        xindex::Int64,yindex::Int64,tindex::Int64,
        nx::Int64,ny::Int64,nt::Int64,Us::U,
        x::T, y::T, t::T,
        )::T where {T,U}
    xp::SVector{4,T} = SVector{4,T}((1.0,x,x^2,x^3))
    yp::SVector{4,T} = SVector{4,T}((1.0,y,y^2,y^3))
    tp::SVector{4,T} = SVector{4,T}((1.0,t,t^2,t^3))
    result::T = zero(T)


    #Specify point values
    @inbounds for k in 0:1, j in 0:1, i in 0:1
        current_index::Tuple{Int64,Int64,Int64} = (xindex + i,yindex + j,tindex + k)
        res::T = Us[earthIndex(current_index,nx,ny,nt)]
        aindex::Int64 = i + 2*j + 4*k +1
        tmpresult::T = zero(T)
        @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
            rowIndex::Int64 = A.rowval[r] - 1
            xValIndex::Int64 = rowIndex % 4
            yValIndex::Int64 = (div( (rowIndex - xValIndex),4)) % 4
            tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)
            tmpresult += A.nzval[r]*xp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
        end
        result += res*tmpresult
    end

    #Specify first derivatives
    @inbounds for k in 0:1, j in 0:1, i in 0:1
        for ddirection in 1:3
            if ddirection == 1
                ddirectionT::Tuple{Int64,Int64,Int64} = (1,0,0)
            elseif ddirection == 2
                ddirectionT = (0,1,0)
            elseif ddirection == 3
                ddirectionT = (0,0,1)
            end
            current_index = (xindex + i,yindex+ j,tindex + k)
            res::T = 0.0
            res += Us[earthIndex(current_index .+ ddirectionT,nx,ny,nt)]
            res -= Us[earthIndex(current_index .- ddirectionT,nx,ny,nt)]
            res /= 2.0
            aindex::Int64 = i + 2*j + 4*k + 8*ddirection + 1

            tmpresult::T = zero(T)

            @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
                rowIndex::Int64 = A.rowval[r] - 1
                xValIndex::Int64 = rowIndex % 4
                yValIndex::Int64 = ( div((rowIndex - xValIndex),4)) % 4
                tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)
                tmpresult += A.nzval[r]*xp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
            end
            result += tmpresult*res
        end
    end

    #Specify (mixed) second derivatives

    @inbounds for k in 0:1, j in 0:1, i in 0:1
        for ddirection1 in 1:2
            if ddirection1 == 1
                ddirection1T::Tuple{Int64,Int64,Int64} = (1,0,0)
            else
                ddirection1T = (0,1,0)
            end
            for ddirection2 in (ddirection1+1):3
                if ddirection2 == 2
                    ddirection2T::Tuple{Int64,Int64,Int64} = (0,1,0)
                else
                    ddirection2T = (0,0,1)
                end
                    current_index::Tuple{Int64,Int64,Int64} = (xindex + i,yindex+j,tindex+k)
                    res = 0.0
                    res += Us[earthIndex(current_index .+ ddirection1T .+ ddirection2T,nx,ny,nt)]
                    res -= Us[earthIndex(current_index .+ ddirection1T .- ddirection2T,nx,ny,nt)]
                    res -= Us[earthIndex(current_index .- ddirection1T .+ ddirection2T,nx,ny,nt)]
                    res += Us[earthIndex(current_index .- ddirection1T .- ddirection2T,nx,ny,nt)]
                    res /= 4.0
                    aindex = i + 2*j + 4*k + (2*(ddirection1-1) + (ddirection2 - ddirection1 - 1) + 4)*8 + 1

                    tmpresult = zero(T)

                    @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
                        rowIndex::Int64 = A.rowval[r] - 1
                        xValIndex::Int64 = rowIndex % 4
                        yValIndex::Int64 = ( div((rowIndex - xValIndex),4)) % 4
                        tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)
                        tmpresult += A.nzval[r]*xp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
                    end
                    result += res*tmpresult
            end
        end
    end
    #Specfiy (mixed) third derivatives
    @inbounds for j in 0:1, k in 0:1, i in 0:1
        current_index = (xindex+i,yindex+j,tindex+k)
        res = 0.0
        res += Us[earthIndex(current_index .+ (1,1,1) ,nx,ny,nt)]
        res -= Us[earthIndex(current_index .+ (1,1,-1) ,ny,ny,nt)]
        res -= Us[earthIndex(current_index .+ (1,-1,1),nx,ny,nt )]
        res += Us[earthIndex(current_index .+ (1,-1,-1),nx,ny,nt )]
        res -= Us[earthIndex(current_index .+ (-1,1,1),nx,ny,nt )]
        res += Us[earthIndex(current_index .+ (-1,1,-1),nx,ny,nt )]
        res += Us[earthIndex(current_index .+ (-1,-1,1),nx,ny,nt )]
        res -= Us[earthIndex(current_index .+ (-1,-1,-1),nx,ny,nt )]
        res /= 8.0
        aindex = i + 2*j + 4*k + 56 + 1

        tmpresult = zero(T)

        @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
            rowIndex = A.rowval[r] - 1
            xValIndex::Int64 = rowIndex % 4
            yValIndex::Int64 = ( div((rowIndex - xValIndex),4)) % 4
            tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)
            tmpresult += A.nzval[r]*xp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
        end
        result += res*tmpresult
    end
    return result
end

@inbounds @inline function base_tricubic_interpolation_gradient(
        xindex::Int64,yindex::Int64,tindex::Int64,
        nx::Int64,ny::Int64,nt::Int64,Us::U,
        x::T, y::T, t::T,
        )::SVector{2,T} where {T,U}
    xp::SVector{4,T} = SVector{4,T}((1.0,x,x^2,x^3))
    dxp::SVector{4,T} = SVector{4,T}((0.0,1.0,2*x,3*x^2))
    yp::SVector{4,T} = SVector{4,T}((1.0,y,y^2,y^3))
    dyp::SVector{4,T} = SVector{4,T}((0.0,1.0,2*y,3*y^2))
    tp::SVector{4,T} = SVector{4,T}((1.0,t,t^2,t^3))
    result1::T = zero(T)
    result2::T = zero(T)

    #Specify point values
    @inbounds for k in 0:1, j in 0:1, i in 0:1
        current_index::Tuple{Int64,Int64,Int64} = (xindex + i,yindex + j,tindex + k)
        res::T = Us[earthIndex(current_index,nx,ny,nt)]
        aindex::Int64 = i + 2*j + 4*k +1
        tmpresult1::T = zero(T)
        tmpresult2::T = zero(T)
        @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
            rowIndex::Int64 = A.rowval[r] - 1
            xValIndex::Int64 = rowIndex % 4
            yValIndex::Int64 = (div( (rowIndex - xValIndex),4)) % 4
            tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)
            tmpresult1 += A.nzval[r]*dxp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
            tmpresult2 += A.nzval[r]*xp[xValIndex + 1]*dyp[yValIndex + 1]*tp[tValIndex + 1]
        end
        result1 += res*tmpresult1
        result2 += res*tmpresult2
    end

    #Specify first derivatives
    @inbounds for k in 0:1, j in 0:1, i in 0:1
        for ddirection in 1:3
            if ddirection == 1
                ddirectionT::Tuple{Int64,Int64,Int64} = (1,0,0)
            elseif ddirection == 2
                ddirectionT = (0,1,0)
            elseif ddirection == 3
                ddirectionT = (0,0,1)
            end
            current_index = (xindex + i,yindex+ j,tindex + k)
            res::T = 0.0
            res += Us[earthIndex(current_index .+ ddirectionT,nx,ny,nt)]
            res -= Us[earthIndex(current_index .- ddirectionT,nx,ny,nt)]
            res /= 2.0
            aindex::Int64 = i + 2*j + 4*k + 8*ddirection + 1

            tmpresult1::T = zero(T)
            tmpresult2::T = zero(T)

            @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
                rowIndex::Int64 = A.rowval[r] - 1
                xValIndex::Int64 = rowIndex % 4
                yValIndex::Int64 = ( div((rowIndex - xValIndex),4)) % 4
                tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)
                tmpresult1 += A.nzval[r]*dxp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
                tmpresult2 += A.nzval[r]*xp[xValIndex + 1]*dyp[yValIndex + 1]*tp[tValIndex + 1]
            end

            result1 += res*tmpresult1
            result2 += res*tmpresult2
        end
    end


    #Specify (mixed) second derivatives

    @inbounds for k in 0:1, j in 0:1, i in 0:1
        @simd for ddirection1 in 1:2
            if ddirection1 == 1
                ddirection1T::Tuple{Int64,Int64,Int64} = (1,0,0)
            else
                ddirection1T = (0,1,0)
            end
            for ddirection2 in (ddirection1+1):3
                if ddirection2 == 2
                    ddirection2T::Tuple{Int64,Int64,Int64} = (0,1,0)
                else
                    ddirection2T = (0,0,1)
                end
                    current_index::Tuple{Int64,Int64,Int64} = (xindex + i,yindex+j,tindex+k)
                    res = 0.0
                    res += Us[earthIndex(current_index .+ ddirection1T .+ ddirection2T,nx,ny,nt)]
                    res -= Us[earthIndex(current_index .+ ddirection1T .- ddirection2T,nx,ny,nt)]
                    res -= Us[earthIndex(current_index .- ddirection1T .+ ddirection2T,nx,ny,nt)]
                    res += Us[earthIndex(current_index .- ddirection1T .- ddirection2T,nx,ny,nt)]
                    res /= 4.0
                    aindex = i + 2*j + 4*k + (2*(ddirection1-1) + (ddirection2 - ddirection1 - 1) + 4)*8 + 1

                    tmpresult1::T = zero(T)
                    tmpresult2::T = zero(T)

                    @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
                        rowIndex::Int64 = A.rowval[r] - 1
                        xValIndex::Int64 = rowIndex % 4
                        yValIndex::Int64 = ( div((rowIndex - xValIndex),4)) % 4
                        tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)

                        tmpresult1 += A.nzval[r]*dxp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
                        tmpresult2 += A.nzval[r]*xp[xValIndex + 1]*dyp[yValIndex + 1]*tp[tValIndex + 1]
                    end

                    result1 += res*tmpresult1
                    result2 += res*tmpresult2
            end
        end
    end
    #Specfiy (mixed) third derivatives
    @inbounds for j in 0:1, k in 0:1, i in 0:1
        current_index = (xindex+i,yindex+j,tindex+k)
        res = 0.0
        res += Us[earthIndex(current_index .+ (1,1,1) ,nx,ny,nt)]
        res -= Us[earthIndex(current_index .+ (1,1,-1) ,ny,ny,nt)]
        res -= Us[earthIndex(current_index .+ (1,-1,1),nx,ny,nt )]
        res += Us[earthIndex(current_index .+ (1,-1,-1),nx,ny,nt )]
        res -= Us[earthIndex(current_index .+ (-1,1,1),nx,ny,nt )]
        res += Us[earthIndex(current_index .+ (-1,1,-1),nx,ny,nt )]
        res += Us[earthIndex(current_index .+ (-1,-1,1),nx,ny,nt )]
        res -= Us[earthIndex(current_index .+ (-1,-1,-1),nx,ny,nt )]
        res /= 8.0
        aindex = i + 2*j + 4*k + 56 + 1

        tmpresult1::T = zero(T)
        tmpresult2::T = zero(T)

        @simd for r in  A.colptr[aindex]:(A.colptr[aindex + 1] - 1)
            rowIndex = A.rowval[r] - 1
            xValIndex::Int64 = rowIndex % 4
            yValIndex::Int64 = ( div((rowIndex - xValIndex),4)) % 4
            tValIndex::Int64 = div(( div((rowIndex - xValIndex),4) - yValIndex) , 4)

            tmpresult1 += A.nzval[r]*dxp[xValIndex + 1]*yp[yValIndex + 1]*tp[tValIndex + 1]
            tmpresult2 += A.nzval[r]*xp[xValIndex + 1]*dyp[yValIndex + 1]*tp[tValIndex + 1]
        end

        result1 += res*tmpresult1
        result2 += res*tmpresult2

    end
    return SVector{2,T}((result1,result2))
end

#Function barrier version
function fast_tricubic_earth_interpolate(u::StaticVector{2,T},p,tin::Float64) where {T<:Real}
    Us = p[1]
    Vs = p[2]
    return fast_tricubic_earth_interpolate_internal(u,Us,Vs,p,tin)
end

#Actual interpolation version
function fast_tricubic_earth_interpolate_internal(u::StaticVector{2,T},Us::S,Vs::U,p,tin::Float64) where {T<:Real,S,U}
    nx::Int64 = size(Us)[1]
    ny::Int64 = size(Us)[2]
    nt::Int64 = size(Us)[3]
    #Get the spatial bounds from p
    ll1::Float64,ur1::Float64  = p[3]
    ll2::Float64,ur2::Float64 = p[4]
    t0::Float64,tf::Float64 = p[5]
    t::Float64 = tin
    if t > tf
        t = tf
    end
    if t < t0
        t = t0
    end
    xindex::Int64, xcoord::T = gooddivrem((mod((u[1] - ll1), 360)*nx)/360.0,1)
    yindex::Int64, ycoord::T = gooddivrem((mod((u[2] - ll2), 180)*ny)/180.0,1)
    tindex::Int64, tcoord::T = gooddivrem((nt-1)*(t-t0)/(tf-t0),1)

    res1::T = base_tricubic_interpolation(xindex,yindex,tindex,nx,ny,nt,Us, xcoord, ycoord,tcoord)
    res2::T = base_tricubic_interpolation(xindex,yindex,tindex,nx,ny,nt,Vs, xcoord, ycoord,tcoord)

    return SVector{2,T}((res1,res2))
end


function fast_tricubic_earth_interpolate_gradient(u::StaticVector{2,T},p,tin::Float64) where {T<:Real}
    sshs = p[7]
    return fast_tricubic_earth_interpolate_gradient_internal(u,sshs,p,tin)
end


function fast_tricubic_earth_interpolate_gradient_internal(u::StaticVector{2,T},sshs::S,p,tin::Float64) where {T<:Real,S,U}
    nx::Int64 = size(sshs)[1]
    ny::Int64 = size(sshs)[2]
    nt::Int64 = size(sshs)[3]
    #Get the spatial bounds from p
    ll1::Float64,ur1::Float64  = p[3]
    ll2::Float64,ur2::Float64 = p[4]
    t0::Float64,tf::Float64 = p[5]
    t::Float64 = tin
    if t > tf
        t = tf
    end
    if t < t0
        t = t0
    end
    xindex::Int64, xcoord::T = gooddivrem((mod((u[1] - ll1), 360)*nx)/360.0,1)
    yindex::Int64, ycoord::T = gooddivrem((mod((u[2] - ll2), 180)*ny)/180.0,1)
    tindex::Int64, tcoord::T = gooddivrem((nt-1)*(t-t0)/(tf-t0),1)

    res1::SVector{2,T} = base_tricubic_interpolation_gradient(xindex,yindex,tindex,nx,ny,nt,sshs, xcoord, ycoord,tcoord)
    return SVector{2,T}((res1[1]*nx/360,res1[2]*ny/180))
end

function interp_inbounds(x,y,p)
    return !isnan(p[3][x,y])
end
function inbounds_checker_bilinear(x,y,p)
    return !isnan(fast_trilinear_earth_interpolate(SVector{2,Float64}((x,y)),p[6],p[6][5][1])[1])
end

function sshVelocityRHS(u,p,t)
    du::SVector{2,Float64} = fast_tricubic_earth_interpolate_gradient(u,p,t)

    g::Float64 = 9.807 #Gravitational constant (m/s^2)
    R::Float64 = 6371e3 #Radius of earth (m)
    Ω = 7.2921159e-5 #Mean angular velocity (rad/s)


    #Rescale to m/radian
    du1 = 360/(2π)*du[1]
    du2 = 360/(2π)*du[2]

    #Calculate velocity
    C = g/(R^2*2*Ω*sin(deg2rad(u[2]))*cos(deg2rad(u[2])))
    #Go from radians/s to radians/day
    C *= 24*3600

    #Go from  radians/day to deg/day
    C *= 360/(2π)

    return  SVector{2,Float64}((-du2*C,du1*C))
end




#=
function fast_trilinear_ssh_gradient_flipped(du::AbstractVector{Float64},u::AbstractVector{Float64},p,tin)
    sshs = p[7]
    nx::Int64 = size(sshs)[1]
    ny::Int64 = size(sshs)[2]
    nt::Int64 = size(sshs)[3]
    #Get the spatial bounds from p
    ll1::Float64,ur1::Float64  = p[3]
    ll2::Float64,ur2::Float64 = p[4]
    t0::Float64,tf::Float64 = p[5]
    t::Float64 = tin
    if t > tf
        t = tf
    end
    if t < t0
        t = t0
    end
    #Just divrem, but casts the first result to Int
    function gooddivrem(x,y)
        a,b = divrem(x,y)
        return Int(a), b
    end

    xindex::Int64, xcoord::Float64 = gooddivrem((mod((u[1] - ll1), 360)*nx)/360.0,1)
    yindex::Int64, ycoord::Float64 = gooddivrem((mod((u[2] - ll2), 180)*ny)/180.0,1)
    #
    tindex::Int64, tcoord::Float64 = gooddivrem((nt-1)*(t-t0)/(tf-t0),1)
    tindex += 1
    #Make sure we don't go out of bounds
    tpp = tindex + 1
    if tpp > nt
        tpp = nt
    end
    #Actual interpolation for u
    r1::Float64 =  sshs[xindex+1,yindex + 1,tindex ]*(1 - xcoord) +            sshs[ (xindex + 1) % nx + 1,yindex + 1, tindex ]*xcoord
    dr1dx::Float64 =  -sshs[xindex+1,yindex + 1,tindex ] +            sshs[ (xindex + 1) % nx + 1,yindex + 1, tindex ]
    r2::Float64 =  sshs[xindex+1,(yindex + 1) % ny + 1,tindex ]*(1 - xcoord) + sshs[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tindex ]*xcoord
    dr2dx::Float64 =  -sshs[xindex+1,(yindex + 1) % ny + 1,tindex ] + sshs[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tindex ]
    r3::Float64 =  sshs[xindex+1,yindex + 1,tpp ]*(1 - xcoord) +               sshs[ (xindex + 1) % nx + 1,yindex + 1, tpp ]*xcoord
    dr3dx::Float64 =  -sshs[xindex+1,yindex + 1,tpp ] +               sshs[ (xindex + 1) % nx + 1,yindex + 1, tpp ]
    r4::Float64 =  sshs[xindex+1,(yindex + 1) % ny + 1,tpp ]*(1 - xcoord) +    sshs[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tpp ]*xcoord
    dr4dx::Float64 =  -sshs[xindex+1,(yindex + 1) % ny + 1,tpp ] +    sshs[ (xindex + 1) % nx + 1,(yindex + 1)%ny + 1, tpp ]

    res1::T = ((1-tcoord)*((1-ycoord)*dr1dx + ycoord*dr2dx)
             + tcoord*((1-ycoord)*dr3dx + ycoord*dr4dx))*nx/360.0

    res2::T = ((1-tcoord)*(-r1+ r2)
             + tcoord*(-r3+ r4))*ny/180.0
     return SVector{2,T}((res1,res2))
end






function intSSHRHS(du,u,p,t)
    sshI = p[1]
    Interpolations.gradient!(du,sshI, u[1],u[2],t)
    tmp::Float64 = du[1]
    du[2] = du[1]
    du[1] = tmp

    g::Float64 = 9.807 #Gravitational constant (m/s^2)
    R::Float64 = 6371e3 #Radius of earth (m)
    Ω = 7.2921159e-5 #Mean angular velocity (rad/s)


    #Rescale to m/radian
    du[2] *= 360/(2π)
    du[1] *= 360/(2π)

    #Calculate velocity
    C = g/(R^2*2*Ω*sin(deg2rad(u[2]))*cos(deg2rad(u[2])))
    #Go from radians/s to radians/day
    C *= 24*3600

    #Go from  radians/day to deg/day
    C *= 360/(2π)

    du[1] *= -C
    du[2] *= C

    return
end

function interpolateSSHPeriodic(Lon,Lat,Time, sshs, interpolation_type)
    # convert arrays into linspace-form for interpolation
    lon = range(minimum(Lon),stop=maximum(Lon),length=length(Lon))
    lat = range(minimum(Lat),stop=maximum(Lat),length=length(Lat))
    time = range(minimum(Time),stop=maximum(Time),length=length(Time))
    sshI = Interpolations.interpolate(sshs,interpolation_type,OnGrid())
    sshI = Interpolations.scale(sshI,lon,lat,time)
    sshE = extrapolate(sshI,(Periodic(),Periodic(),Flat()))
    return sshE
end
=#
