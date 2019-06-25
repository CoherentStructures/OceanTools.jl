
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


"""
    build_full_finite_difference_matrix()

Returns a (sparse) matrix representation of the matrix that calculates
the vector `b` from Lekien & Marsden's paper using first-order finite differences
from grid values.

"""
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

    for ddirection in 1:3
        for i in 0:1, j in 0:1, k in 0:1
            if ddirection == 1
                ddirectionT::Tuple{Int64,Int64,Int64} = (1,0,0)
            elseif ddirection == 2
                ddirectionT = (0,1,0)
            elseif ddirection == 3
                ddirectionT = (0,0,1)
            end
            aindex::Int64 = i + 2*j + 4*k + 8*ddirection + 1
            result[aindex, to1d(i + ddirectionT[1], j + ddirectionT[2], k+ddirectionT[3])] = 4
            result[aindex, to1d(i - ddirectionT[1], j - ddirectionT[2], k-ddirectionT[3])] = -4
        end
    end

    #Specify (mixed) second derivatives
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

                    aindex = i + 2*j + 4*k + 32 +
                        8*((ddirection1-1) + (ddirection2 - ddirection1 - 1)) + 1
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
    getIndex(x,x0,xf,nx,boundary_behaviour)

Calculates the indexes `i,j` and local coordinates corresponding to a real number `x`
where `x` is in c_i, c_j and the interval [x0,xf) is partitioned into intervals
[x0 = c_0, c_1), [c_1,c_2), ... [c_(nx-1), xf) where each intervals has equal length.
If `boundary_behaviour==0`, periodic boundary behaviour is assumed.
If `boundary_behaviour==1`, flat boundary behaviour is assumed.
If `boundary_behaviour==2`, an out of bounds error is thrown is x is outside of the interval.
"""
@inline function getIndex(
    x::Float64, x0::Float64,xf::Float64,
    nx::Int64,boundary_behaviour::Int
    )::Tuple{Int64,Int64,Float64}
    if boundary_behaviour == 0 #Periodic boundary
        xindex::Int64, xcoord::Float64 = gooddivrem((mod(x - x0, (xf-x0))*(nx))/(xf-x0),1)
        xpp = (xindex+1) % nx
    elseif boundary_behaviour >= 1 # Flat boundary (==1) or Error (==2)
        xindex, xcoord = gooddivrem(((x-x0)*nx)/(xf-x0),1)
        xpp = xindex + 1
        if xpp >= nx
            if boundary_behaviour==2
                throw(BoundsError("Out of bounds access"))
            else
                xpp = (nx-1)
            end
        end
        if xindex >= nx
            xindex = (nx-1)
        end
        if xindex < 0
            if boundary_behaviour==2
                throw(BoundsError("Out of bounds access"))
            else
                xindex = 0
            end
        end
        if xpp < 0
            xpp = 0
        end
    else
        throw(AssertionError("Bad value for boundary_behaviour"))
    end
    return xindex, xpp, xcoord
end


@inline function getIndex2(
    x::Float64, x0::Float64,xf::Float64,
    nx::Int64,boundary_behaviour::Int
    )::Tuple{Int64,Int64,Int64,Int64,Float64}
    if boundary_behaviour == 0 #Periodic boundary
        xindex::Int64, xcoord::Float64 = gooddivrem((mod(x - x0, (xf-x0))*(nx))/(xf-x0),1)
        xpp = (xindex+1) % nx
        xpp2 = (xindex+2) % nx
        xmm = mod((xindex-1),nx)
    elseif boundary_behaviour >= 1 # Flat boundary (==1) or Error (==2)
        xindex, xcoord = gooddivrem(((x-x0)*nx)/(xf-x0),1)
        xpp = xindex + 1
        xpp2 = xindex + 2
        xmm = xindex-1

        if xpp2 >= nx
            if boundary_behaviour==2
                throw(BoundsError("Out of bounds access"))
            else
                xpp2 = (nx-1)
            end
        end
        if xpp >= nx
            xpp = (nx-1)
        end
        if xindex >= nx
            xindex = (nx-1)
        end
        if xmm >= nx
            xmm = (nx-1)
        end

        if xmm < 0
            if boundary_behaviour==2
                throw(BoundsError("Out of bounds access"))
            else
                xmm = 0
            end
        end

        if xindex < 0
                xindex = 0
        end
        if xpp < 0
            xpp = 0
        end

        if xpp2 < 0
            xpp2 = 0
        end

    else
        throw(AssertionError("Bad value for boundary_behaviour"))
    end
    return xindex, xpp,xpp2,xmm, xcoord
end

struct ItpMetadata{T}
    nx::Int
    ny::Int
    nt::Int
    LL::SVector{3,Float64}
    UR::SVector{3,Float64}
    data::T
    boundaryX::Int64
    boundaryY::Int64
    boundaryT::Int64

    function ItpMetadata(nx::Int,ny::Int,nt::Int,
                        LL::AbstractArray{Float64},UR::AbstractArray{Float64},
                        data::T, boundaryX::Int64,
                        boundaryY::Int64,boundaryT::Int64
                        ) where T
        @assert length(LL) == 3
        @assert length(UR) == 3
        new{T}(nx,ny,nt, (@SVector [LL[1],LL[2],LL[3]]),
                        (@SVector [UR[1],UR[2],UR[3]]),
                        data,boundaryX,boundaryY,boundaryT
                        )
    end
end


"""
    uv_trilinear(u,p,tin)

Trilinear interpolation of velocity field at `u` at time `tin`.
Velocity field stored in `p` as returned by `getP`
Periodic boundary in x and y, constant in t direction
"""
function uv_trilinear(
        u::SArray{Tuple{2},T,1,2},p::ItpMetadata{S},tin::Float64
        )::SArray{Tuple{2},T,1,2} where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return uv_trilinear_internal(u,p,tin,Us,Vs)
end


function uv_trilinear_internal(
        u::SArray{Tuple{2},T,1,2},p::ItpMetadata{S},tin::Float64,Us::U,Vs::U
        )::SArray{Tuple{2},T,1,2} where {T<:Real,S,U}

    #Get data from p
    nx::Int64, ny::Int64, nt::Int64 = p.nx,p.ny,p.nt
    ll1::Float64,ll2::Float64,t0::Float64 = p.LL
    ur1::Float64,ur2::Float64,tf::Float64 = p.UR

    @inbounds xindex::Int64,xpp::Int64, xcoord::T = getIndex(u[1],ll1,ur1, nx, p.boundaryX)
    @inbounds yindex::Int64,ypp::Int64, ycoord::T = getIndex(u[2],ll2,ur2, ny, p.boundaryY)
    tindex::Int64,tpp::Int64, tcoord::T = getIndex(tin,t0,tf, nt, p.boundaryT)

    @inbounds begin
    r1u::T =  Us[xindex+1,yindex + 1,tindex + 1 ]*(1 - xcoord) +
                    Us[ xpp + 1,yindex + 1, tindex + 1 ]*xcoord
    r2u::T =  Us[xindex+1,ypp + 1,tindex + 1 ]*(1 - xcoord) +
                    Us[ xpp + 1,ypp + 1, tindex + 1 ]*xcoord
    r3u::T =  Us[xindex + 1,yindex + 1,tpp + 1 ]*(1 - xcoord) +
                    Us[ xpp + 1,yindex + 1, tpp + 1 ]*xcoord
    r4u::T =  Us[xindex+1,ypp + 1,tpp + 1 ]*(1 - xcoord) +
                    Us[ xpp + 1,ypp + 1, tpp + 1 ]*xcoord
    res1::T =  (
        (1-tcoord)*((1-ycoord)*r1u + ycoord*r2u)
         + tcoord*((1-ycoord)*r3u + ycoord*r4u))

    r1v::T =  Vs[xindex+1,yindex + 1,tindex + 1 ]*(1 - xcoord) +
                    Vs[ xpp + 1,yindex + 1, tindex + 1 ]*xcoord
    r2v::T =  Vs[xindex+1,ypp + 1,tindex + 1 ]*(1 - xcoord) +
                    Vs[ xpp + 1,ypp + 1, tindex + 1 ]*xcoord
    r3v::T =  Vs[xindex + 1,yindex + 1,tpp + 1 ]*(1 - xcoord) +
                    Vs[ xpp + 1,yindex + 1, tpp + 1 ]*xcoord
    r4v::T =  Vs[xindex+1,ypp + 1,tpp + 1 ]*(1 - xcoord) +
                    Vs[ xpp + 1,ypp + 1, tpp + 1 ]*xcoord
    res2::T =  (
        (1-tcoord)*((1-ycoord)*r1v + ycoord*r2v)
         + tcoord*((1-ycoord)*r3v + ycoord*r4v))
     end

    return SArray{Tuple{2},T,1,2}((res1,res2))
end

@inline function base_tricubic_interpolation(
        xindex::Int64,yindex::Int64,tindex::Int64,
        xpp::Int64,ypp::Int64,tpp::Int64,
        xpp2::Int64,ypp2::Int64,tpp2::Int64,
        xmm::Int64,ymm::Int64,tmm::Int64,
        nx::Int64,ny::Int64,
        x::T, y::T, t::T,
        Us::U,
        Vs::U
        )::SVector{2,T} where {T,U}
    xp::SVector{4,T} = SVector{4,T}((1.0,x,x^2,x^3))
    yp::SVector{4,T} = SVector{4,T}((1.0,y,y^2,y^3))
    tp::SVector{4,T} = SVector{4,T}((1.0,t,t^2,t^3))

    function earthIndexRaw(i::Int64,j::Int64,k::Int64)::Int64
        xi::Int64 = 0
        if i == -1
            xi = xmm
        elseif i == 0
            xi = xindex
        elseif i == 1
            xi = xpp
        elseif i == 2
            xi = xpp2
        end
        yi::Int64 = 0

        if j == -1
            yi = ymm
        elseif j == 0
            yi = yindex
        elseif j == 1
            yi = ypp
        elseif j == 2
            yi = ypp2
        end

        ti::Int64 = 0
        if k == -1
            ti = tmm
        elseif k == 0
            ti = tindex
        elseif k == 1
            ti = tpp
        elseif k == 2
            ti = tpp2
        end
        return xi + yi * nx + ti*nx*ny + 1
    end

    uvals = @inbounds @SArray T[Us[earthIndexRaw(i,j,k)] for i in -1:2, j in -1:2, k in -1:2]
    @inbounds AFbyu = A_times_svec(F_times_svec(uvals))
    vvals = @inbounds @SArray T[Vs[earthIndexRaw(i,j,k)] for i in -1:2, j in -1:2, k in -1:2]
    @inbounds AFbyv = A_times_svec(F_times_svec(vvals))
    #@inbounds AFbyv = A_times_svec(
    #SVector{64,Float64}(build_full_finite_difference_matrix()*Vector{Float64}(vec(vvals)))
    #)

    res1::T = zero(T)
    res2::T = zero(T)
    @inbounds for i in 1:4,j in 1:4, k in 1:4
            res1 += xp[i]*yp[j]*tp[k]*AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
            res2 += xp[i]*yp[j]*tp[k]*AFbyv[(i-1) + 4*(j-1) + 16*(k-1) + 1]
    end
    return @SVector [res1/8,res2/8]
end

@inline function base_tricubic_interpolation_gradient(
        xindex::Int64,yindex::Int64,tindex::Int64,
        xpp::Int64,ypp::Int64,tpp::Int64,
        xpp2::Int64,ypp2::Int64,tpp2::Int64,
        xmm::Int64,ymm::Int64,tmm::Int64,
        nx::Int64,ny::Int64,
        x::T, y::T, t::T,
        Us::U,
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
        xi::Int64 = 0
        if i == -1
            xi = xmm
        elseif i == 0
            xi = xindex
        elseif i == 1
            xi = xpp
        elseif i == 2
            xi = xpp2
        end
        yi::Int64 = 0

        if j == -1
            yi = ymm
        elseif j == 0
            yi = yindex
        elseif j == 1
            yi = ypp
        elseif j == 2
            yi = ypp2
        end

        ti::Int64 = 0
        if k == -1
            ti = tmm
        elseif k == 0
            ti = tindex
        elseif k == 1
            ti = tpp
        elseif k == 2
            ti = tpp2
        end
        return xi + yi * nx + ti*nx*ny + 1
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
function uv_tricubic(u::StaticVector{2,T},p::ItpMetadata{S},tin::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return uv_tricubic_internal(u,Us,Vs,p,tin)
end

#Actual interpolation version
function uv_tricubic_internal(u::StaticVector{2,T},Us::U,Vs::U,p::ItpMetadata{S},tin::Float64) where {T<:Real,S,U}

    nx::Int64, ny::Int64, nt::Int64 = p.nx,p.ny,p.nt
    ll1::Float64,ll2::Float64,t0::Float64 = p.LL
    ur1::Float64,ur2::Float64,tf::Float64 = p.UR

    @inbounds xindex::Int64,xpp::Int64,xpp2::Int64,xmm::Int64, xcoord::T = getIndex2(u[1],ll1,ur1, nx, p.boundaryX)
    @inbounds yindex::Int64,ypp::Int64,ypp2::Int64,ymm::Int64, ycoord::T = getIndex2(u[2],ll2,ur2, ny, p.boundaryY)
    tindex::Int64,tpp::Int64,tpp2::Int64,tmm::Int64, tcoord::T = getIndex2(tin,t0,tf, nt, p.boundaryT)

    return base_tricubic_interpolation(
        xindex,yindex,tindex,
        xpp,ypp,tpp,
        xpp2,ypp2,tpp2,
        xmm,ymm,tmm,
        nx,ny,
        xcoord, ycoord,tcoord,
        Us,Vs
        )
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

    nx::Int64, ny::Int64, nt::Int64 = p.nx,p.ny,p.nt
    ll1::Float64,ll2::Float64,t0::Float64 = p.LL
    ur1::Float64,ur2::Float64,tf::Float64 = p.UR

    @inbounds xindex::Int64,xpp::Int64,xpp2::Int64,xmm::Int64, xcoord::T = getIndex2(u[1],ll1,ur1, nx, p.boundaryX)
    @inbounds yindex::Int64,ypp::Int64,ypp2::Int64,ymm::Int64, ycoord::T = getIndex2(u[2],ll2,ur2, ny, p.boundaryY)
    tindex::Int64,tpp::Int64,tpp2::Int64,tmm::Int64, tcoord::T = getIndex2(tin,t0,tf, nt, p.boundaryT)

    return base_tricubic_interpolation(
        xindex,yindex,tindex,
        xpp,ypp,tpp,
        xpp2,ypp2,tpp2,
        xmm,ymm,tmm,
        nx,ny,
        xcoord, ycoord,tcoord,
        sshs,sshs
        )[1] #TODO: Modify this to avoid doing the work twice

end



function ssh_tricubic_gradient(u::StaticVector{2,T},p::ItpMetadata{S},tin::Float64) where {T<:Real,S}
    sshs = p.data
    return ssh_tricubic_gradient_internal(u,sshs,p,tin)
end


function ssh_tricubic_gradient_internal(u::StaticVector{2,T},sshs::U,p::ItpMetadata{S},tin::Float64) where {T<:Real,S,U}

    nx::Int64, ny::Int64, nt::Int64 = p.nx,p.ny,p.nt
    ll1::Float64,ll2::Float64,t0::Float64 = p.LL
    ur1::Float64,ur2::Float64,tf::Float64 = p.UR

    @inbounds xindex::Int64,xpp::Int64,xpp2::Int64,xmm::Int64, xcoord::T = getIndex2(u[1],ll1,ur1, nx, p.boundaryX)
    @inbounds yindex::Int64,ypp::Int64,ypp2::Int64,ymm::Int64, ycoord::T = getIndex2(u[2],ll2,ur2, ny, p.boundaryY)
    tindex::Int64,tpp::Int64,tpp2::Int64,tmm::Int64, tcoord::T = getIndex2(tin,t0,tf, nt, p.boundaryT)

    res1 = base_tricubic_interpolation_gradient(
        xindex,yindex,tindex,
        xpp,ypp,tpp,
        xpp2,ypp2,tpp2,
        xmm,ymm,tmm,
        nx,ny,
        xcoord, ycoord,tcoord,
        sshs
        )

    return SVector{2,T}((res1[1]*nx/(ur1-ll1),res1[2]*ny/(ur2-ll2)))
end


#Function barrier version
function uv_tricubic_eqvari(u::StaticMatrix{2,3,T},p::ItpMetadata{S},tin::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return uv_tricubic_eqvari_internal(u,Us,Vs,p,tin)
end


function uv_tricubic_eqvari_internal(
    uIn::StaticMatrix{2,3,T}, Us::U,Vs::U,p::ItpMetadata{S},tin::Float64) where {T<:Real,S,U}
    u = @SVector [uIn[1,1], uIn[2,1]]


    nx::Int64, ny::Int64, nt::Int64 = p.nx,p.ny,p.nt
    ll1::Float64,ll2::Float64,t0::Float64 = p.LL
    ur1::Float64,ur2::Float64,tf::Float64 = p.UR

    @inbounds xindex::Int64,xpp::Int64,xpp2::Int64,xmm::Int64, xcoord::T = getIndex2(u[1],ll1,ur1, nx, p.boundaryX)
    @inbounds yindex::Int64,ypp::Int64,ypp2::Int64,ymm::Int64, ycoord::T = getIndex2(u[2],ll2,ur2, ny, p.boundaryY)
    tindex::Int64,tpp::Int64,tpp2::Int64,tmm::Int64, tcoord::T = getIndex2(tin,t0,tf, nt, p.boundaryT)

    Uitp::SVector{3,T} = base_tricubic_interpolation_gradient(
        xindex,yindex,tindex,
        xpp,ypp,tpp,
        xpp2,ypp2,tpp2,
        xmm,ymm,tmm,
        nx,ny,
        xcoord, ycoord,tcoord,
        Us
        )

    Vitp::SVector{3,T} = base_tricubic_interpolation_gradient(
        xindex,yindex,tindex,
        xpp,ypp,tpp,
        xpp2,ypp2,tpp2,
        xmm,ymm,tmm,
        nx,ny,
        xcoord, ycoord,tcoord,
        Vs
        )

    return @SMatrix [Uitp[3] (Uitp[1]*uIn[1,2]*nx/(ur1-ll1) + Uitp[2]*uIn[2,2]*ny/(ur2-ll2)) (Uitp[1]*uIn[1,3]*nx/(ur1-ll1) + Uitp[2]*uIn[2,3]*ny/(ur2-ll2));
    Vitp[3] (Vitp[1]*uIn[1,2]*nx/(ur1-ll1)  + Vitp[2]*uIn[2,2]*ny/(ur2-ll2)) (Vitp[1]*uIn[1,3]*nx/(ur1-ll1) + Vitp[2]*uIn[2,3]*ny/(ur2-ll2))]
end

function interp_inbounds(x,y,p)
    return !isnan(p[3][x,y])
end
function inbounds_checker_bilinear(x,y,p)
    return !isnan(uv_trilinear(SVector{2,Float64}((x,y)),p[6],p[6][5][1])[1])
end

function sshVelocityRHS(u,p::ItpMetadata{S},t::Float64) where S
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
