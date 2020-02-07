
#Constants needed for tricubic interpolation
include("coeff.jl")

@enum BoundaryBehaviour periodic=0 flat=1 outofbounds=2


function sparse_by_svec(A::SparseMatrixCSC{TA}, x::Symbol) where {TA}
    n, m = size(A)

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

@generated function A_times_svec(x::SArray)
    sparse_by_svec(getA(), :x)
end

@generated function F_times_svec(x::SArray)
    sparse_by_svec(build_full_finite_difference_matrix(), :x)
end

@generated function AF_times_svec(x::SArray)
    tobuild = getA() * build_full_finite_difference_matrix()
    sparse_by_svec(tobuild, :x)
end


"""
    build_full_finite_difference_matrix()

Return a (sparse) matrix representation of the matrix that calculates the vector `b` from
Lekien & Marsden's paper using first-order finite differences from grid values.
"""
function build_full_finite_difference_matrix()
    result = zeros(Int, 64, 64)
    function to1d(i, j, k)
        return (i+1) + 4*(j+1) + 16*(k+1) + 1
    end

    #Point values
    for i in 0:1, j in 0:1, k in 0:1
        aindex = i + 2j + 4k + 1
        result[aindex, to1d(i, j, k)] = 8
    end

    #Specify first derivatives

    for ddirection in 1:3
        for i in 0:1, j in 0:1, k in 0:1
            if ddirection == 1
                ddirectionT = (1, 0, 0)
            elseif ddirection == 2
                ddirectionT = (0, 1, 0)
            elseif ddirection == 3
                ddirectionT = (0, 0, 1)
            end
            aindex = i + 2j + 4k + 8ddirection + 1
            result[aindex, to1d(i+ddirectionT[1], j+ddirectionT[2], k+ddirectionT[3])] =  4
            result[aindex, to1d(i-ddirectionT[1], j-ddirectionT[2], k-ddirectionT[3])] = -4
        end
    end

    #Specify (mixed) second derivatives

    for i in 0:1, j in 0:1, k in 0:1
        for ddirection1 in 1:2
            if ddirection1 == 1
                ddirection1T = (0, 1, 0)
            else
                ddirection1T = (1, 0, 0)
            end
            for ddirection2 in (ddirection1+1):3
                if ddirection2 == 2
                    ddirection2T = (1, 0, 0)
                else
                    ddirection2T = (0, 0, 1)
                end
                    shift = 0
                    if ddirection1 == 1
                        if ddirection2 == 2
                            #xy derivative
                            shift = 0
                        else
                            #xz derivative
                            shift = 2
                        end
                    else
                        #yz derivative
                        shift = 1
                    end

                    aindex = i + 2j + 4k + 32 + 8shift + 1
                    result[aindex, to1d(
                             i + ddirection1T[1] + ddirection2T[1],
                             j + ddirection1T[2] + ddirection2T[2],
                             k + ddirection1T[3] + ddirection2T[3]
                             )] = 2

                    result[aindex, to1d(
                            i + ddirection1T[1] - ddirection2T[1],
                            j + ddirection1T[2] - ddirection2T[2],
                            k + ddirection1T[3] - ddirection2T[3]
                            )] = -2

                    result[aindex, to1d(
                            i - ddirection1T[1] + ddirection2T[1],
                            j - ddirection1T[2] + ddirection2T[2],
                            k - ddirection1T[3] + ddirection2T[3]
                            )] = -2

                    result[aindex, to1d(
                            i - ddirection1T[1] - ddirection2T[1],
                            j - ddirection1T[2] - ddirection2T[2],
                            k - ddirection1T[3] - ddirection2T[3]
                            )] = 2
            end
        end
    end

    #Specfiy (mixed) third derivatives
    for i in 0:1, j in 0:1, k in 0:1
        aindex = i + 2j + 4k + 56 + 1

        result[aindex, to1d(i+1,j+1,k+1)] =  1
        result[aindex, to1d(i+1,j+1,k-1)] = -1
        result[aindex, to1d(i+1,j-1,k+1)] = -1
        result[aindex, to1d(i+1,j-1,k-1)] =  1
        result[aindex, to1d(i-1,j+1,k+1)] = -1
        result[aindex, to1d(i-1,j+1,k-1)] =  1
        result[aindex, to1d(i-1,j-1,k+1)] =  1
        result[aindex, to1d(i-1,j-1,k-1)] = -1
    end
    return sparse(result)
end

const A = getA()
const F = build_full_finite_difference_matrix()
const AF = A*build_full_finite_difference_matrix()



#Just divrem, but casts the first result to Int
@inline function gooddivrem(x::T, y::Int64)::Tuple{Int64,T} where {T}
    a, b = divrem(x, T(y))
    return Base.unsafe_trunc(Int,a), b
end

#=
function gooddivrem(x::ForwardDiff.Dual, y)
    a,b = divrem(x,y)
    return Int(ForwardDiff.value(a)), b
end
=#

"""
    getIndex(x, x0, xf, nx, boundary_behaviour)

Calculates the indexes `i,j` and local coordinates corresponding to a real number `x`
where `x` is in c_i, c_j and the interval [x0,xf) is partitioned into intervals
[x0 = c_0, c_1), [c_1,c_2), ... [c_(nx-1), xf) where all intervals have equal length.
"""
@inline function getIndex(x::Float64, x0::Float64, xf::Float64, nx::Int64,
                        boundary_behaviour::BoundaryBehaviour)::Tuple{Int64,Int64,Float64}
    if boundary_behaviour == periodic #Periodic boundary
        xindex, xcoord = gooddivrem((mod(x - x0, (xf-x0))*(nx))/(xf-x0), 1)
        xpp = (xindex+1) % nx
    elseif boundary_behaviour == flat || boundary_behaviour == outofbounds
        xindex, xcoord = gooddivrem(((x-x0)*nx)/(xf-x0), 1)
        xpp = xindex + 1
        if xpp >= nx
            if boundary_behaviour == outofbounds
                throw(BoundsError("Out of bounds access"))
            else
                xpp = (nx-1)
            end
        end
        if xindex >= nx
            xindex = (nx-1)
        end
        if xindex < 0
            if boundary_behaviour == outofbounds
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


@inline function getIndex2(x::Float64, x0::Float64, xf::Float64, nx::Int64,
                boundary_behaviour::BoundaryBehaviour)::Tuple{Int64,Int64,Int64,Int64,Float64}
    if boundary_behaviour == periodic #Periodic boundary
        xindex, xcoord = gooddivrem((mod(x - x0, (xf-x0))*(nx))/(xf-x0), 1)
        xindex = xindex % nx
        xpp = (xindex+1) % nx
        xpp2 = (xindex+2) % nx
        xmm = mod((xindex-1), nx)
    elseif boundary_behaviour == flat || boundary_behaviour == outofbounds # Flat boundary (==1) or Error (==2)
        xindex, xcoord = gooddivrem(((x-x0)*nx)/(xf-x0), 1)
        xpp = xindex + 1
        xpp2 = xindex + 2
        xmm = xindex-1

        if xpp2 >= nx
            if boundary_behaviour == outofbounds
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
            if boundary_behaviour == outofbounds
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
    return xindex, xpp, xpp2, xmm, xcoord
end


struct ItpMetadata{T}
    nx::Int
    ny::Int
    nt::Int
    LL::SVector{3,Float64}
    UR::SVector{3,Float64}
    data::T
    boundaryX::BoundaryBehaviour
    boundaryY::BoundaryBehaviour
    boundaryT::BoundaryBehaviour

    function ItpMetadata(nx::Int, ny::Int, nt::Int,
                        LL::AbstractArray{Float64}, UR::AbstractArray{Float64}, data::T,
                        boundaryX::BoundaryBehaviour,
                        boundaryY::BoundaryBehaviour,
                        boundaryT::BoundaryBehaviour) where {T}
        @assert length(LL) == 3
        @assert length(UR) == 3
        new{T}(nx, ny, nt, (@SVector [LL[1], LL[2], LL[3]]), (@SVector [UR[1],UR[2],UR[3]]),
                        data, boundaryX, boundaryY, boundaryT)
    end

    function ItpMetadata(nx::Int, ny::Int, nt::Int,
                        LL::AbstractArray{Float64}, UR::AbstractArray{Float64}, data::T,
                        boundaryX::Int64,
                        boundaryY::Int64,
                        boundaryT::Int64
                        ) where T
        @assert length(LL) == 3
        @assert length(UR) == 3
        new{T}(nx, ny, nt, (@SVector [LL[1],LL[2],LL[3]]), (@SVector [UR[1],UR[2],UR[3]]), data,
            BoundaryBehaviour(boundaryX), BoundaryBehaviour(boundaryY), BoundaryBehaviour(boundaryT))
    end

    function ItpMetadata(
                        xspan::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
                        yspan::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
                        tspan::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
                        data::T,
                        boundaryX::BoundaryBehaviour,
                        boundaryY::BoundaryBehaviour,
                        boundaryT::BoundaryBehaviour) where {T}
        nx = length(xspan)
        ny = length(yspan)
        nt = length(tspan)
        LL = (@SVector [minimum(xspan),minimum(yspan),minimum(tspan)])
        UR = (@SVector [maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(tspan)+step(tspan)])
        new{T}(nx, ny, nt, LL,UR, data, boundaryX, boundaryY, boundaryT)
    end

end

struct ItpMetadata3{T}
    nx::Int
    ny::Int
    nz::Int
    nt::Int
    LL::SVector{4,Float64}
    UR::SVector{4,Float64}
    data::T
    boundaryX::Int
    boundaryY::Int
    boundaryZ::Int
    boundaryT::Int

    function ItpMetadata3(nx::Int, ny::Int, nz::Int, nt::Int,
                        LL::AbstractArray{Float64},UR::AbstractArray{Float64}, data::T,
                        boundaryX::Int, boundaryY::Int,boundaryZ::Int,boundaryT::Int) where {T}
        @assert length(LL) == 4
        @assert length(UR) == 4
        new{T}(nx, ny, nz, nt, (@SVector [LL[1], LL[2], LL[3], LL[4]]),
                        (@SVector [UR[1], UR[2], UR[3], UR[3]]), data,
                        boundaryX, boundaryY, boundaryZ, boundaryT)
    end
end



"""
    uv_trilinear(u, p, t)

Trilinear interpolation of velocity field at `u` at time `t`.
Velocity field stored in `p` as returned by [`getP`](@ref).
Periodic boundary in x and y, constant in t direction.
"""
function uv_trilinear(u::SVector{2,T}, p::ItpMetadata{S}, t::Float64)::SArray{Tuple{2},T,1,2} where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return _uv_trilinear(u, Us, Vs, p, t)
end


@inbounds function _uv_trilinear(u::SVector{2,T}, Us::U, Vs::U, p::ItpMetadata{S}, t::Float64) where {T<:Real,S,U}
    #Get data from p
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR

    xindex, xpp, xcoord = getIndex(u[1], ll1, ur1, nx, p.boundaryX)
    yindex, ypp, ycoord = getIndex(u[2], ll2, ur2, ny, p.boundaryY)
    tindex, tpp, tcoord = getIndex(t, t0, tf, nt, p.boundaryT)

    r1u = Us[xindex + 1, yindex + 1, tindex + 1]*(1-xcoord) + Us[xpp + 1, yindex + 1, tindex + 1]*xcoord
    r2u = Us[xindex + 1, ypp + 1, tindex + 1]*(1-xcoord) + Us[xpp + 1, ypp + 1, tindex + 1]*xcoord
    r3u = Us[xindex + 1, yindex + 1, tpp + 1]*(1-xcoord) + Us[xpp + 1, yindex + 1, tpp + 1]*xcoord
    r4u = Us[xindex + 1, ypp + 1, tpp + 1]*(1-xcoord) + Us[xpp + 1,ypp + 1, tpp + 1 ]*xcoord
    res1 = ((1-tcoord)*((1-ycoord)*r1u + ycoord*r2u) + tcoord*((1-ycoord)*r3u + ycoord*r4u))

    r1v = Vs[xindex+1, yindex + 1, tindex + 1]*(1 - xcoord) + Vs[xpp + 1, yindex + 1, tindex + 1]*xcoord
    r2v = Vs[xindex+1, ypp + 1, tindex + 1]*(1 - xcoord) + Vs[xpp + 1, ypp + 1, tindex + 1]*xcoord
    r3v = Vs[xindex + 1, yindex + 1,tpp + 1]*(1 - xcoord) + Vs[xpp + 1, yindex + 1, tpp + 1]*xcoord
    r4v = Vs[xindex+1, ypp + 1, tpp + 1]*(1 - xcoord) + Vs[xpp + 1, ypp + 1, tpp + 1]*xcoord
    res2 = ((1-tcoord)*((1-ycoord)*r1v + ycoord*r2v) + tcoord*((1-ycoord)*r3v + ycoord*r4v))

    return SVector{2,T}((res1, res2))
end

@inline function base_tricubic_interpolation(
        xindex::Int, yindex::Int, tindex::Int,
        xpp::Int, ypp::Int, tpp::Int,
        xpp2::Int, ypp2::Int, tpp2::Int,
        xmm::Int, ymm::Int, tmm::Int,
        nx::Int, ny::Int,
        x::T, y::T, t::T,
        Us::U, Vs::U) where {T,U}
    xp = SVector{4,T}((1.0, x, x^2, x^3))
    yp = SVector{4,T}((1.0, y, y^2, y^3))
    tp = SVector{4,T}((1.0, t, t^2, t^3))

    function earthIndexRaw(i::Int, j::Int, k::Int)
        if i == -1
            xi = xmm
        elseif i == 0
            xi = xindex
        elseif i == 1
            xi = xpp
        elseif i == 2
            xi = xpp2
        else
            xi = 0
        end

        if j == -1
            yi = ymm
        elseif j == 0
            yi = yindex
        elseif j == 1
            yi = ypp
        elseif j == 2
            yi = ypp2
        else
            yi = 0
        end

        if k == -1
            ti = tmm
        elseif k == 0
            ti = tindex
        elseif k == 1
            ti = tpp
        elseif k == 2
            ti = tpp2
        else
            ti = 0
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

    res1 = zero(T)
    res2 = zero(T)
    @inbounds for i in 1:4, j in 1:4, k in 1:4
            res1 += xp[i]*yp[j]*tp[k]*AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
            res2 += xp[i]*yp[j]*tp[k]*AFbyv[(i-1) + 4*(j-1) + 16*(k-1) + 1]
    end
    return SVector{2,T}((res1/8.0, res2/8.0))
end

@inline function base_tricubic_interpolation_gradient(
        xindex::Int, yindex::Int, tindex::Int,
        xpp::Int, ypp::Int, tpp::Int,
        xpp2::Int, ypp2::Int, tpp2::Int,
        xmm::Int, ymm::Int, tmm::Int,
        nx::Int, ny::Int,
        x::T, y::T, t::T,
        Us::U) where {T,U}

    xp = SVector{4,T}((1.0, x, x^2, x^3))
    dxp = SVector{4,T}((0.0, 1.0, 2*x, 3*x^2))
    yp = SVector{4,T}((1.0, y, y^2, y^3))
    dyp = SVector{4,T}((0.0, 1.0, 2*y, 3*y^2))
    tp = SVector{4,T}((1.0, t, t^2, t^3))
    result1 = zero(T)
    result2 = zero(T)
    result3 = zero(T)

    function earthIndexRaw(i::Int, j::Int, k::Int)
        if i == -1
            xi = xmm
        elseif i == 0
            xi = xindex
        elseif i == 1
            xi = xpp
        elseif i == 2
            xi = xpp2
        else
            xi = 0
        end

        if j == -1
            yi = ymm
        elseif j == 0
            yi = yindex
        elseif j == 1
            yi = ypp
        elseif j == 2
            yi = ypp2
        else
            yi = 0
        end

        if k == -1
            ti = tmm
        elseif k == 0
            ti = tindex
        elseif k == 1
            ti = tpp
        elseif k == 2
            ti = tpp2
        else
            ti = 0
        end
        return xi + yi * nx + ti*nx*ny + 1
    end

    uvals = @inbounds @SArray T[Us[earthIndexRaw(i,j,k)] for i in -1:2, j in -1:2, k in -1:2]

    @inbounds AFbyu = A_times_svec(F_times_svec(uvals))
    @inbounds for i in 1:4, j in 1:4, k in 1:4
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
    return SVector{3,T}((result1/8.0, result2/8.0, result3/8.0))
end

"""
    uv_tricubic(x, p, t)

Tricubic interpolation of velocity field at `u` at time `t`.
Velocity field stored in `p` as returned by [`getP`](@ref).
Periodic boundary in x and y, constant in t direction.
"""
function uv_tricubic(u::SVector{2,T}, p::ItpMetadata{S}, t::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return _uv_tricubic(u, Us, Vs, p, t)
end

#Actual interpolation version
function _uv_tricubic(u::SVector{2,T}, Us::U, Vs::U, p::ItpMetadata{S}, t::Float64) where {T<:Real,S,U}

    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR

    @inbounds xindex, xpp, xpp2, xmm, xcoord = getIndex2(u[1], ll1, ur1, nx, p.boundaryX)
    @inbounds yindex, ypp, ypp2, ymm, ycoord = getIndex2(u[2], ll2, ur2, ny, p.boundaryY)
    @inbounds tindex, tpp, tpp2, tmm, tcoord = getIndex2(t, t0, tf, nt, p.boundaryT)

    return base_tricubic_interpolation(
        xindex, yindex, tindex,
        xpp, ypp, tpp,
        xpp2, ypp2, tpp2,
        xmm, ymm, tmm,
        nx, ny,
        xcoord, ycoord, tcoord,
        Us, Vs)
end



"""
    ssh_tricubic(x, p, t)
Tricubic interpolation of scalar field field at `u` at time `t`.
Scalar field stored in `p` as returned by [`getP`](@ref).
Periodic boundary in x and y, constant in t direction.
"""
function ssh_tricubic(u::StaticVector{2,T}, p, t::Float64) where {T<:Real}
    sshs = p[7]
    return _ssh_tricubic(u, sshs, p, t)
end


function _ssh_tricubic(u::StaticVector{2,T}, sshs::S, p, t::Float64) where {T<:Real,S}

    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR

    @inbounds xindex, xpp, xpp2, xmm, xcoord = getIndex2(u[1], ll1, ur1, nx, p.boundaryX)
    @inbounds yindex, ypp, ypp2, ymm, ycoord = getIndex2(u[2], ll2, ur2, ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, tcoord = getIndex2(t, t0, tf, nt, p.boundaryT)

    return base_tricubic_interpolation(
        xindex, yindex, tindex,
        xpp, ypp, tpp,
        xpp2, ypp2, tpp2,
        xmm, ymm, tmm,
        nx, ny,
        xcoord, ycoord, tcoord,
        sshs, sshs)[1] #TODO: Modify this to avoid doing the work twice
end

function ssh_tricubic_gradient(u::StaticVector{2,T}, p::ItpMetadata{S}, t::Float64) where {T<:Real,S}
    sshs = p.data
    return _ssh_tricubic_gradient(u, sshs, p, t)
end


function _ssh_tricubic_gradient(u::StaticVector{2,T}, sshs::U, p::ItpMetadata{S}, t::Float64) where {T<:Real,S,U}

    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR

    @inbounds xindex, xpp, xpp2, xmm, xcoord = getIndex2(u[1],ll1,ur1, nx, p.boundaryX)
    @inbounds yindex, ypp, ypp2, ymm, ycoord = getIndex2(u[2],ll2,ur2, ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, tcoord = getIndex2(t, t0, tf, nt, p.boundaryT)

    res1 = base_tricubic_interpolation_gradient(
        xindex, yindex, tindex,
        xpp, ypp, tpp,
        xpp2, ypp2, tpp2,
        xmm, ymm, tmm,
        nx, ny,
        xcoord, ycoord, tcoord,
        sshs)

    return SVector{2,T}((res1[1]*nx/(ur1-ll1), res1[2]*ny/(ur2-ll2)))
end


#Function barrier version
function uv_tricubic_eqvari(u::StaticMatrix{2,3,T}, p::ItpMetadata{S}, t::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return _uv_tricubic_eqvari(u, Us, Vs, p, t)
end


function _uv_tricubic_eqvari(uIn::StaticMatrix{2,3,T}, Us::U, Vs::U, p::ItpMetadata{S}, t::Float64) where {T<:Real,S,U}
    u = @SVector [uIn[1,1], uIn[2,1]]

    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR

    @inbounds xindex, xpp, xpp2, xmm, xcoord = getIndex2(u[1], ll1, ur1, nx, p.boundaryX)
    @inbounds yindex, ypp, ypp2, ymm, ycoord = getIndex2(u[2], ll2, ur2, ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, tcoord = getIndex2(t, t0, tf, nt, p.boundaryT)

    Uitp = base_tricubic_interpolation_gradient(
        xindex, yindex, tindex,
        xpp, ypp, tpp,
        xpp2, ypp2, tpp2,
        xmm, ymm, tmm,
        nx, ny,
        xcoord, ycoord, tcoord,
        Us)

    Vitp = base_tricubic_interpolation_gradient(
        xindex, yindex, tindex,
        xpp, ypp, tpp,
        xpp2, ypp2, tpp2,
        xmm, ymm, tmm,
        nx, ny,
        xcoord, ycoord, tcoord,
        Vs)

    return @SMatrix [Uitp[3] (Uitp[1]*uIn[1,2]*nx/(ur1-ll1) + Uitp[2]*uIn[2,2]*ny/(ur2-ll2)) (Uitp[1]*uIn[1,3]*nx/(ur1-ll1) + Uitp[2]*uIn[2,3]*ny/(ur2-ll2));
    Vitp[3] (Vitp[1]*uIn[1,2]*nx/(ur1-ll1)  + Vitp[2]*uIn[2,2]*ny/(ur2-ll2)) (Vitp[1]*uIn[1,3]*nx/(ur1-ll1) + Vitp[2]*uIn[2,3]*ny/(ur2-ll2))]
end

function interp_inbounds(x, y, p)
    return !isnan(p[3][x, y])
end
function inbounds_checker_bilinear(x, y, p)
    return !isnan(uv_trilinear(SVector{2,Float64}((x,y)), p[6], p[6][5][1])[1])
end

function sshVelocityRHS(u, p::ItpMetadata{S}, t::Float64) where S
    ∇h = ssh_tricubic_gradient(u, p, t)

    g = 9.807 #Gravitational constant (m/s^2)
    R = 6371e3 #Radius of earth (m)
    Ω = 7.2921159e-5 #Mean angular velocity (rad/s)

    # convert from m/deg to m/rad
    scale = 360/(2π)
    ∇h *= scale

    # calculate velocity conversion factor in rad/s
    C = g/(R^2*2*Ω*sin(deg2rad(u[2]))*cos(deg2rad(u[2])))
    # convert from rad/s to rad/day
    C *= 24*3600
    # convert from rad/day to deg/day
    C *= scale

    return  SVector{2,Float64}((-∇h[2]*C, ∇h[1]*C))
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

"""
    uv_quadlinear(u, p, t)

Quadlinear interpolation of velocity field at `u` at time `t`.
Velocity field stored in `p` as returned by [`getP`](@®ef).
Periodic boundary in x and y, constant in t direction.
"""
function uv_quadlinear(u::SVector{3,T}, p::ItpMetadata3{S}, t::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    Ws = p.data[3]
    return _uv_quadlinear(u, Us, Vs, Ws, p, t)
end

function _uv_quadlinear(u::SVector{3,T}, Us::U, Vs::U, Ws::U, p::ItpMetadata3{S}, t::Float64) where {T<:Real,S,U}
    #Get data from p
    nx, ny, nz, nt = p.nx, p.ny, p.nz, p.nt
    ll1, ll2, ll3, t0 = p.LL
    ur1, ur2, ur3, tf = p.UR

    @inbounds xindex, xpp, xcoord = getIndex(u[1], ll1, ur1, nx, p.boundaryX)
    @inbounds yindex, ypp, ycoord = getIndex(u[2], ll2, ur2, ny, p.boundaryY)
    @inbounds zindex, zpp, zcoord = getIndex(u[3], ll3, ur3, nz, p.boundaryZ)
    tindex, tpp, tcoord = getIndex(t, t0, tf, nt, p.boundaryT)

    @inbounds begin
        function singleTriLinear(what,zindex)
            r1u = what[xindex+1, yindex + 1, zindex + 1, tindex + 1]*(1 - xcoord) +
                            what[xpp + 1, yindex + 1, zindex + 1, tindex + 1]*xcoord
            r2u = what[xindex+1, ypp + 1, zindex + 1, tindex + 1]*(1 - xcoord) +
                            what[xpp + 1, ypp + 1, zindex + 1, tindex + 1]*xcoord
            r3u = what[xindex + 1, yindex + 1, zindex + 1, tpp + 1]*(1 - xcoord) +
                            what[xpp + 1, yindex + 1, zindex + 1, tpp + 1]*xcoord
            r4u = what[xindex+1, ypp + 1, zindex + 1, tpp + 1]*(1 - xcoord) +
                            what[xpp + 1, ypp + 1, zindex + 1, tpp + 1]*xcoord
            return ((1-tcoord)*((1-ycoord)*r1u + ycoord*r2u) + tcoord*((1-ycoord)*r3u + ycoord*r4u))
        end
        function singleQuadLinear(what)
            return (1-zcoord)*singleTriLinear(what, zindex) + zcoord*singleTriLinear(what, zpp)
        end
    end
    return SVector{3,T}((singleQuadLinear(Us), singleQuadLinear(Vs), singleQuadLinear(Ws)))
end
