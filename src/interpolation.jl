
#Constants needed for tricubic interpolation
include("coeff.jl")

function sparse_by_svec(A::SparseMatrixCSC{TA}, x::Symbol) where TA
    n,m = size(A)

    coeff_initializers = Expr(:block, [:($(Symbol("y",i)) = zero(T)) for i in 1:n]...)
    coeff_updates = Expr(:block)
    rows = rowvals(A)
    vals = nonzeros(A)
    for col in 1:m
        for j in nzrange(A,col)
            row = rows[j]
            val = vals[j]
            push!(coeff_updates.args, :($(Symbol("y",row)) += $val *$x[$col]))
        end
    end
    return quote
        @inbounds begin
            T = promote_type($TA, eltype($x))
            $coeff_initializers
            $coeff_updates
            $(:(return SVector(tuple($([Symbol("y",i) for i in 1:n]...)))))
        end
    end
end

@generated function A_times_svec(x::SArray{m}) where m
    sparse_by_svec(getA(),:x)
end

@generated function F_times_svec(x::SArray{m}) where m
    sparse_by_svec(build_full_finite_difference_matrix(),:x)
end

@generated function AF_times_svec(x::SArray{m}) where m
    tobuild = getA()*build_full_finite_difference_matrix()
    sparse_by_svec(tobuild,:x)
end


function build_full_finite_difference_matrix()
    result = spzeros(Int64,64,64)
    function to1d(i,j,k)
        return (i+1) + 4*(j+1) + 16*(k+1) + 1
    end

    #Point values
    for i in 0:1, j in 0:1, k in 0:1
        aindex::Int64 = i + 2*j + 4*k +1
        result[aindex,to1d(i,j,k)] = 8
    end

    #Specify first derivatives
    for i in 0:1, j in 0:1, k in 0:1
        for ddirection in 1:3
            if ddirection == 1
                ddirectionT::Tuple{Int64,Int64,Int64} = (1,0,0)
            elseif ddirection == 2
                ddirectionT = (0,1,0)
            elseif ddirection == 3
                ddirectionT = (0,0,1)
            end
            #res += Us[earthIndex(current_index .+ ddirectionT,nx,ny,nt)]
            #res -= Us[earthIndex(current_index .- ddirectionT,nx,ny,nt)]
            aindex::Int64 = i + 2*j + 4*k + 8*ddirection + 1
            result[aindex, to1d(i + ddirectionT[1], j + ddirectionT[2], k+ddirectionT[3])] = 4
            result[aindex, to1d(i - ddirectionT[1], j - ddirectionT[2], k-ddirectionT[3])] = -4
        end
    end

    for i in 0:1, j in 0:1, k in 0:1
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

                    aindex = i + 2*j + 4*k + (2*(ddirection1-1) + (ddirection2 - ddirection1 - 1) + 4)*8 + 1
                    result[aindex,to1d(
                             i + ddirection1T[1] + ddirection2T[1],
                             j + ddirection1T[2] + ddirection2T[2],
                             k + ddirection1T[3] + ddirection2T[3]
                             )] = 2

                    result[aindex,to1d(
                            i + ddirection1T[1] - ddirection2T[1],
                            j + ddirection1T[2] - ddirection2T[2],
                            k + ddirection1T[3] - ddirection2T[3]
                            )] = -2

                    result[aindex,to1d(
                            i - ddirection1T[1] + ddirection2T[1],
                            j - ddirection1T[2] + ddirection2T[2],
                            k - ddirection1T[3] + ddirection2T[3]
                            )] = -2

                    result[aindex,to1d(
                            i - ddirection1T[1] - ddirection2T[1],
                            j - ddirection1T[2] - ddirection2T[2],
                            k - ddirection1T[3] - ddirection2T[3]
                            )] = 2
            end
        end
    end
    #Specfiy (mixed) third derivatives
    for i in 0:1, j in 0:1, k in 0:1
        aindex = i + 2*j + 4*k + 56 + 1

        res = 0.0
        result[aindex,to1d(i+1,j+1,k+1)] = 1
        result[aindex,to1d(i+1,j+1,k-1)] = -1
        result[aindex,to1d(i+1,j-1,k+1)] = -1
        result[aindex,to1d(i+1,j-1,k-1)] = 1
        result[aindex,to1d(i-1,j+1,k+1)] = -1
        result[aindex,to1d(i-1,j+1,k-1)] = 1
        result[aindex,to1d(i-1,j-1,k+1)] = 1
        result[aindex,to1d(i-1,j-1,k-1)] = -1
    end
    return result
end

const A = getA()
const FT = permutedims(build_full_finite_difference_matrix())
const F = build_full_finite_difference_matrix()
const AF = A*build_full_finite_difference_matrix()



#Just divrem, but casts the first result to Int
@inline function gooddivrem(x::T,y::Int64)::Tuple{Int64,T} where T
    a,b = divrem(x,T(y))
    return Base.unsafe_trunc(Int,a), b
end

#=
function gooddivrem(x::ForwardDiff.Dual, y)
    a,b = divrem(x,y)
    return Int(ForwardDiff.value(a)), b
end
=#

"""
    uv_trilinear(u,p,tin)

Trilinear interpolation of velocity field at `u` at time `tin`.
Velocity field stored in `p` as returned by `getP`
Periodic boundary in x and y, constant in t direction
"""
function uv_trilinear(
        u::SArray{Tuple{2},T,1,2},p,tin::Float64
        )::SArray{Tuple{2},T,1,2} where {T<:Real}
    Us = p[1]
    Vs = p[2]
    return uv_trilinear_internal(u,p,tin,Us,Vs)
end

function uv_trilinear_internal(
        u::SArray{Tuple{2},T,1,2},p,tin::Float64,Us::U,Vs::V
        )::SArray{Tuple{2},T,1,2} where {T<:Real,U,V}

    nx::Int64, ny::Int64, nt::Int64 = size(Us)
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
    #TODO: Maybe rewrite this with the earthIndex function (see below)
    @inbounds begin
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
    end
    return SArray{Tuple{2},T,1,2}((res1,res2))
end

#Periodic boundary in x and y directions, constant boundary in z direction.
@inline @inbounds function earthIndex(x::Tuple{Int64,Int64,Int64},nx::Int64,ny::Int64,nt::Int64)::CartesianIndex{3}
    return CartesianIndex( (mod(x[1],nx) + 1 ,mod(x[2],ny) + 1, max(min(x[3], nt-1 ),0) + 1))
end

#Some basic tricubic interpolation
@inline function base_tricubic_interpolation(
        xindex::Int64,yindex::Int64,tindex::Int64,
        nx::Int64,ny::Int64,nt::Int64,Us::U,
        x::T, y::T, t::T,
        )::T where {T,U}
    xp::SVector{4,T} = SVector{4,T}((1.0,x,x^2,x^3))
    yp::SVector{4,T} = SVector{4,T}((1.0,y,y^2,y^3))
    tp::SVector{4,T} = SVector{4,T}((1.0,t,t^2,t^3))

    function earthIndexRaw(i::Int64,j::Int64,k::Int64)::Int64
        i_new::Int64 = mod(xindex + i,nx)
        j_new::Int64 = mod(yindex + j,ny)
        k_new::Int64 =  max(min(tindex + k, nt-1 ),0)
        return i_new + j_new * nx + k_new*nx*ny + 1
    end

    uvals = @inbounds @SArray T[Us[earthIndexRaw(i,j,k)] for i in -1:2, j in -1:2, k in -1:2]
    res::T = zero(T)
    @inbounds AFbyu = A_times_svec(F_times_svec(uvals))
    @inbounds for i in 1:4,j in 1:4, k in 1:4
            res += xp[i]*yp[j]*tp[k]*AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
    end
    return res/8.0
end

@inbounds @inline function base_tricubic_interpolation_gradient(
        xindex::Int64,yindex::Int64,tindex::Int64,
        nx::Int64,ny::Int64,nt::Int64,Us::U,
        x::T, y::T, t::T,
        )::SVector{3,T} where {T,U}
    xp::SVector{4,T} = SVector{4,T}((1.0,x,x^2,x^3))
    dxp::SVector{4,T} = SVector{4,T}((0.0,1.0,2*x,3*x^2))
    yp::SVector{4,T} = SVector{4,T}((1.0,y,y^2,y^3))
    dyp::SVector{4,T} = SVector{4,T}((0.0,1.0,2*y,3*y^2))
    tp::SVector{4,T} = SVector{4,T}((1.0,t,t^2,t^3))
    result1::T = zero(T)
    result2::T = zero(T)
    result3::T = zero(T)


    function earthIndexRaw(i::Int64,j::Int64,k::Int64)::Int64
        i_new::Int64 = mod(xindex + i,nx)
        j_new::Int64 = mod(yindex + j,ny)
        k_new::Int64 =  max(min(tindex + k, nt-1 ),0)
        return i_new + j_new * nx + k_new*nx*ny + 1
    end

    uvals = @inbounds @SArray T[Us[earthIndexRaw(i,j,k)] for i in -1:2, j in -1:2, k in -1:2]

    @inbounds AFbyu = A_times_svec(F_times_svec(uvals))
    @inbounds for i in 1:4,j in 1:4, k in 1:4
            curx = xp[i]
            curdx = dxp[i]
            cury = yp[j]
            curdy = dyp[j]
            curt = tp[k]
            aval = AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
            result1 += curdx*cury*curt*aval
            result2 += curx*curdy*curt*aval
            result3 += curx*cury*curt*aval
    end
    return SVector{3,T}((result1/8.0,result2/8.0,result3/8.0))
end

"""
    uv_tricubic(x,p,tin)
Tricubic interpolation of velocity field at `u` at time `tin`.
Velocity field stored in `p` as returned by `getP`
Periodic boundary in x and y, constant in t direction
"""
function uv_tricubic(u::StaticVector{2,T},p,tin::Float64) where {T<:Real}
    Us = p[1]
    Vs = p[2]
    return uv_tricubic_internal(u,Us,Vs,p,tin)
end

#Actual interpolation version
function uv_tricubic_internal(u::StaticVector{2,T},Us::S,Vs::U,p,tin::Float64) where {T<:Real,S,U}
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



"""
    ssh_tricubic(x,p,tin)
Tricubic interpolation of scalar field field at `u` at time `tin`.
Scalar field stored in `p` as returned by `getP`
Periodic boundary in x and y, constant in t direction
"""
function ssh_tricubic(u::StaticVector{2,T},p,tin::Float64) where {T<:Real}
    sshs = p[7]
    return ssh_tricubic_internal(u,sshs,p,tin)
end


function ssh_tricubic_internal(u::StaticVector{2,T},sshs::S,p,tin::Float64) where {T<:Real,S}
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

    res1::T = base_tricubic_interpolation(xindex,yindex,tindex,nx,ny,nt,sshs, xcoord, ycoord,tcoord)

    return res1
end



function ssh_tricubic_gradient(u::StaticVector{2,T},p,tin::Float64) where {T<:Real}
    sshs = p[7]
    return ssh_tricubic_gradient_internal(u,sshs,p,tin)
end


function ssh_tricubic_gradient_internal(u::StaticVector{2,T},sshs::S,p,tin::Float64) where {T<:Real,S,U}
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

    res1::SVector{3,T} = base_tricubic_interpolation_gradient(xindex,yindex,tindex,nx,ny,nt,sshs, xcoord, ycoord,tcoord)
    return SVector{2,T}((res1[1]*nx/360,res1[2]*ny/180))
end


#Function barrier version
function uv_tricubic_eqvari(u::StaticMatrix{2,3,T},p,tin::Float64) where {T<:Real}
    Us = p[1]
    Vs = p[2]
    return uv_tricubic_eqvari_internal(u,Us,Vs,p,tin)
end


function uv_tricubic_eqvari_internal(
    uIn::StaticMatrix{2,3,T}, Us::S,Vs::S,p,tin::Float64) where {T<:Real,S,U}
    u = @SVector [uIn[1,1], uIn[2,1]]
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

    Uitp::SVector{3,T} = base_tricubic_interpolation_gradient(xindex,yindex,tindex,nx,ny,nt,Us, xcoord, ycoord,tcoord)
    Vitp::SVector{3,T} = base_tricubic_interpolation_gradient(xindex,yindex,tindex,nx,ny,nt,Vs, xcoord, ycoord,tcoord)

    return @SMatrix [Uitp[3] (Uitp[1]*uIn[1,2]*nx/360 + Uitp[2]*uIn[2,2]*ny/180) (Uitp[1]*uIn[1,3]*nx/360 + Uitp[2]*uIn[2,3]*ny/180);
    Vitp[3] (Vitp[1]*uIn[1,2]*nx/360  + Vitp[2]*uIn[2,2]*ny/180) (Vitp[1]*uIn[1,3]*nx/360 + Vitp[2]*uIn[2,3]*ny/180)]
end

function interp_inbounds(x,y,p)
    return !isnan(p[3][x,y])
end
function inbounds_checker_bilinear(x,y,p)
    return !isnan(uv_trilinear(SVector{2,Float64}((x,y)),p[6],p[6][5][1])[1])
end

function sshVelocityRHS(u,p,t)
    du::SVector{2,Float64} = ssh_tricubic_gradient(u,p,t)

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
