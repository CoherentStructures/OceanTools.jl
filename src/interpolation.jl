
#Constants needed for tricubic interpolation
include("coeff.jl")

@enum BoundaryBehaviour periodic=0 flat=1 outofbounds=2 semiperiodic=3

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
            push!(coeff_updates.args, :($(Symbol("y", row)) += $val *$x[$col]))
        end
    end
    return quote
        @inbounds begin
            T = promote_type($TA, eltype($x))
            $coeff_initializers
            $coeff_updates
            $(:(return SVector(tuple($([Symbol("y", i) for i in 1:n]...)))))
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
@inline function gooddivrem(x::T, y) where {T}
    a, b = divrem(x, convert(T, y))
    return Base.unsafe_trunc(Int, a), b
end

#=
function gooddivrem(x::ForwardDiff.Dual, y)
    a,b = divrem(x,y)
    return Int(ForwardDiff.value(a)), b
end
=#

"""
    getIndex(x, x0, xf, xper, nx, boundary)

Calculates the indexes `i` and `j` and local coordinates corresponding to a real
number `x`, where `x` is in c_i, c_j and the interval [x0,xf) is partitioned into
intervals `[x0 = c_0, c_1), [c_1,c_2), ... [c_(nx-1), xf)` where all intervals
have equal length.
"""
@inline function getIndex(x, x0, xf, xper, nx, boundary::BoundaryBehaviour)
    if boundary == periodic
        #spacing = (xf-x0)/nx
        #xindex, xcoord = gooddivrem(mod(x - x0, xf-x0),spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x - x0, xf - x0) * nx) / (xf - x0), 1)
        xpp = (xindex+1) % nx
    elseif boundary == flat
        #spacing = (xf-x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x - x0)*nx) / (xf - x0), 1)
        xpp = max(min(xindex + 1, nx - 1), 0)
        xindex = max(min(xindex, nx - 1), 0)
    elseif boundary  == outofbounds
        #spacing = (xf - x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x - x0)*nx) / (xf - x0), 1)
        xpp = xindex + 1
        if xpp >= nx || xindex < 0
            error("Out of bounds access at coordinate $x")
        end
    else #boundary == semiperiodic
        #spacing = mod(xf-x0,xper)/nx
        #xindex, xcoord = gooddivrem(mod(x-x0,xper), spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x - x0, xper)*nx) / mod(xf - x0, xper), 1)
        xpp = xindex + 1
        if xpp >= nx || xindex < 0
            error("Out of bounds access at coordinate $x")
        end
    end
    return xindex, xpp, xcoord
end

"""
    getIndex2(args)

Like `getIndex`, but also returns one more set of points to the right/left.
"""
@inline function getIndex2(x, x0, xf, xper, nx, boundary::BoundaryBehaviour)
    if boundary == periodic
        #spacing = (xf-x0)/nx
        #xindex, xcoord = gooddivrem(mod(x - x0, xf-x0),spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x - x0, (xf - x0))*nx) / (xf - x0), 1)
        xindex = xindex % nx
        xpp = (xindex+1) % nx
        xpp2 = (xindex+2) % nx
        xmm = mod(xindex - 1, nx)
    elseif boundary == flat
        #spacing = (xf - x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x - x0)*nx) / (xf - x0), 1)
        xindex = max(min(xindex, nx - 1), 0)
        xpp = max(min(xindex + 1, nx - 1), 0)
        xpp2 = max(min(xindex + 2, nx - 1), 0)
        xmm = max(min(xindex - 1, nx - 1), 0)
    elseif boundary == outofbounds
        #spacing = (xf - x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x - x0)*nx) / (xf - x0), 1)
        xpp = xindex + 1
        xpp2 = xindex + 2
        xmm = xindex - 1
        if xmm < 0 || xpp2 > (nx-1)
            error("Out of bounds access at coordinate $x")
        end
    else #boundary == semiperiodic
        #spacing = mod(xf-x0,xper)/nx
        #xindex, xcoord = gooddivrem(mod(x-x0,xper), spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x - x0, xper)*nx) / mod(xf - x0, xper), 1)
        xpp = xindex + 1
        xpp2 = xindex + 2
        xmm = xindex - 1
        if xmm < 0 || xpp2 > (nx-1)
            error("Out of bounds access at coordinate $x")
        end
    end
    return xmm, xindex, xpp, xpp2, xcoord
end


struct ItpMetadata{T}
    #The number of datapoints in each direction
    nx::Int
    ny::Int
    nt::Int
    #The coordinates of the lower left/upper-right corner in space-time. Note that the upper-right corner must be
    #*one grid spacing beyond* the rightmost datapoint!
    LL::SVector{3,Float64}
    UR::SVector{3,Float64}
    #In the case that a boundary is `semiperiodic`, what the period in this direction is
    periods::SVector{3,Float64}
    #Actual Data
    data::T
    #How to deal with points outside the boundary in each direction
    boundaryX::BoundaryBehaviour
    boundaryY::BoundaryBehaviour
    boundaryT::BoundaryBehaviour
    function ItpMetadata(nx::Int, ny::Int, nt::Int,
                        LL::AbstractVector{Float64}, UR::AbstractVector{Float64}, periods::AbstractVector{Float64}, data::T,
                        boundaryX::BoundaryBehaviour, boundaryY::BoundaryBehaviour, boundaryT::BoundaryBehaviour) where {T}
        @assert length(LL) == 3
        @assert length(UR) == 3
        @assert length(periods) == 3
        if boundaryX != semiperiodic
            @assert periods[1] == 0.0
            @assert UR[1] > LL[1]
        else
            @assert periods[1] != 0
        end

        if boundaryY != semiperiodic
            @assert periods[2] == 0.0
            @assert UR[2] > LL[2]
        else
            @assert periods[2] != 0
        end

        if boundaryT != semiperiodic
            @assert periods[3] == 0.0
            @assert UR[3] > LL[3]
        else
            @assert periods[3] != 0
        end

        return new{T}(
            nx,
            ny,
            nt,
            SVector{3}((LL[1], LL[2], LL[3])),
            SVector{3}((UR[1], UR[2], UR[3])),
            SVector{3}((periods[1],periods[2],periods[3])),
            data,
            boundaryX,
            boundaryY,
            boundaryT,
        )
    end

end


@deprecate ItpMetadata(nx::Int, ny::Int, nt::Int,
                    LL::AbstractVector{Float64}, UR::AbstractVector{Float64}, data::T,
                    boundaryX::Int, boundaryY::Int, boundaryT::Int
                    ) where {T} ItpMetadata(nx, ny, nt, LL, UR,SVector{3}((0.0,0.0,0.0)), data, BoundaryBehaviour(boundaryX), BoundaryBehaviour(boundaryY), BoundaryBehaviour(boundaryT))

function ItpMetadata(xspan::AbstractRange, yspan::AbstractRange, tspan::AbstractRange, data::T,
        boundaryX::BoundaryBehaviour, boundaryY::BoundaryBehaviour, boundaryT::BoundaryBehaviour; periods=SVector{3}((0.0,0.0,0.0))) where {T}
    nx = length(xspan)
    ny = length(yspan)
    nt = length(tspan)
    LL = SVector{3}((minimum(xspan), minimum(yspan), minimum(tspan)))
    UR = SVector{3}((maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(tspan)+step(tspan)))
    return ItpMetadata(nx, ny, nt, LL, UR, periods, data, boundaryX, boundaryY, boundaryT)
end

#For 4D interpolation
struct ItpMetadata3{T}
    nx::Int
    ny::Int
    nz::Int
    nt::Int
    LL::SVector{4,Float64}
    UR::SVector{4,Float64}
    periods::SVector{4,Float64}
    data::T
    boundaryX::BoundaryBehaviour
    boundaryY::BoundaryBehaviour
    boundaryZ::BoundaryBehaviour
    boundaryT::BoundaryBehaviour
end

function ItpMetadata3(nx::Int, ny::Int, nz::Int, nt::Int,
        LL::AbstractVector{Float64}, UR::AbstractVector{Float64}, data::T,
        boundaryX::BoundaryBehaviour, boundaryY::BoundaryBehaviour, boundaryZ::BoundaryBehaviour, boundaryT::BoundaryBehaviour) where {T}
    @assert length(LL) == 4
    @assert length(UR) == 4
    #TODO: implement this!
    @assert boundaryX != semiperiodic
    @assert boundaryY != semiperiodic
    @assert boundaryZ != semiperiodic

    ItpMetadata3(
        nx,
        ny,
        nz,
        nt,
        SVector{4}((LL[1], LL[2], LL[3], LL[4])),
        SVector{4}((UR[1], UR[2], UR[3], UR[4])),
        SVector{4}((0.0,0.0,0.0,0.0)),
        data,
        boundaryX,
        boundaryY,
        boundaryZ,
        boundaryT,
    )
end
@deprecate ItpMetadata3(nx::Int, ny::Int, nz::Int, nt::Int,
        LL::AbstractVector{Float64}, UR::AbstractVector{Float64}, data::T,
        boundaryX::Int, boundaryY::Int, boundaryZ::Int, boundaryT::Int
        ) where {T} ItpMetadata3(nx, ny, nz, nt, LL, UR, data,
                BoundaryBehaviour(boundaryX), BoundaryBehaviour(boundaryY), BoundaryBehaviour(boundaryZ), BoundaryBehaviour(boundaryT))
function ItpMetadata3(xspan::AbstractRange, yspan::AbstractRange, zspan::AbstractRange, tspan::AbstractRange,
        data::T, boundaryX::BoundaryBehaviour, boundaryY::BoundaryBehaviour, boundaryZ::BoundaryBehaviour, boundaryT::BoundaryBehaviour) where {T}
    nx = length(xspan)
    ny = length(yspan)
    nz = length(zspan)
    nt = length(tspan)
    LL = SVector{4}((minimum(xspan), minimum(yspan), minimum(zspan), minimum(tspan)))
    UR = SVector{4}((maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(zspan)+step(zspan), maximum(tspan)+step(tspan)))
    return ItpMetadata3(nx, ny, nz, nt, LL, UR, data, boundaryX, boundaryY, boundaryZ, boundaryT)
end

# TRILINEAR INTERPOLATION

@inline function _get_index(u, p, t)
    #Get data from p
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px, py, pt = p.periods
    x, y = u

    xindex, xpp, xcoord = getIndex(x, ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ycoord = getIndex(y, ll2, ur2, py, ny, p.boundaryY)
    tindex, tpp, tcoord = getIndex(t,  t0,  tf, pt, nt, p.boundaryT)
    
    return xindex, xpp, xcoord, yindex, ypp, ycoord, tindex, tpp, tcoord
end

"""
    uv_trilinear(u, p, t)

Trilinear interpolation of velocity field at `u` at time `t`.
Velocity field stored in `p.data[1]` and `p.data[2]`.
"""
function uv_trilinear(u::SVector{2}, p::ItpMetadata, t)
    inds = _get_index(u, p, t)
    return SVector{2}((_interp_trilinear(p.data[1], inds...), _interp_trilinear(p.data[2], inds...)))
end

"""
    scalar_trilinear(x, p, t)
Trilinear interpolation of scalar field at `x` at time `t`.
Scalar field is assumed to be stored in `p.data[1]`.
"""
function scalar_trilinear(u::SVector{2}, p, t)
    return @inbounds _interp_trilinear(p.data[1], _get_index(u, p, t)...)
end

@inline function _interp_trilinear(Fs, xindex, xpp, xcoord, yindex, ypp, ycoord, tindex, tpp, tcoord)
    r1u = Fs[xindex + 1, yindex + 1, tindex + 1]*(1-xcoord) + Fs[xpp + 1, yindex + 1, tindex + 1]*xcoord
    r2u = Fs[xindex + 1, ypp + 1, tindex + 1]*(1-xcoord) + Fs[xpp + 1, ypp + 1, tindex + 1]*xcoord
    r3u = Fs[xindex + 1, yindex + 1, tpp + 1]*(1-xcoord) + Fs[xpp + 1, yindex + 1, tpp + 1]*xcoord
    r4u = Fs[xindex + 1, ypp + 1, tpp + 1]*(1-xcoord) + Fs[xpp + 1,ypp + 1, tpp + 1 ]*xcoord
    return ((1-tcoord)*((1-ycoord)*r1u + ycoord*r2u) + tcoord*((1-ycoord)*r3u + ycoord*r4u))
end

# TRICUBIC INTERPOLATION

function _get_index2(u::SVector{2}, p, t)
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px, py, pt = p.periods
    xi, yi = u

    xmm, xindex, xpp, xpp2, x = getIndex2(xi, ll1, ur1, px, nx, p.boundaryX)
    ymm, yindex, ypp, ypp2, y = getIndex2(yi, ll2, ur2, py, ny, p.boundaryY)
    tmm, tindex, tpp, tpp2, t = getIndex2(t,   t0,  tf, pt, nt, p.boundaryT)

    xs = (xmm, xindex, xpp, xpp2)
    ys = (ymm, yindex, ypp, ypp2)
    ts = (tmm, tindex, tpp, tpp2)
    
    return xs, x, ys, y, ts, t, nx, ny
end

function _get_cubic_args(u, p, t)
    xs, x, ys, y, ts, t, nx, ny = _get_index2(u, p, t)
    
    xp = SVector{4}((1.0, x, x^2, x^3))
    yp = SVector{4}((1.0, y, y^2, y^3))
    tp = SVector{4}((1.0, t, t^2, t^3))

    return xs, ys, ts, xp, yp, tp, nx, ny
end

"""
    uv_tricubic(u, p, t)

Component-wise tricubic interpolation (Lekien-Marsden + finite differences for values not specified in their paper) of velocity field at `u` at time `t`.
Velocity field stored in `p.data[1]` and `p.data[2]`.
"""
function uv_tricubic(u::SVector{2}, p, t)
    return _interp2_tricubic(p.data[1], p.data[2], _get_cubic_args(u, p, t)...)
end

"""
    scalar_tricubic(x, p, t)

Tricubic interpolation (Lekien-Marsden + finite differences for values not
specified in their paper) of scalar field at `u` at time `t`.
Scalar field is assumed to be stored in `p.data[1]`.
"""
function scalar_tricubic(u::StaticVector{2}, p, t::Float64)
    return _interp1_tricubic(p.data[1], _get_cubic_args(u, p, t)...)
end

@inline _to_raw_index(xi, yj, tk, nx, ny) = xi + nx*(yj + tk*ny) + 1

function _interp1_tricubic(Fs, xs, ys, ts, xp, yp, tp, nx, ny)
    @inbounds begin
        fvals = @SArray [Fs[_to_raw_index(xs[i], ys[j], ts[k], nx, ny)] for i in 1:4, j in 1:4, k in 1:4]
        AFbyu = A_times_svec(F_times_svec(fvals))

        T = eltype(AFbyu)
        res = zero(T)
        for i in 1:4, j in 1:4, k in 1:4
            res += xp[i]*yp[j]*tp[k]*AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
        end
    end
    return res/8.0
end
function _interp2_tricubic(Us, Vs, xs, ys, ts, xp, yp, tp, nx, ny)
    @inbounds begin
        uvals = @SArray [Us[_to_raw_index(xs[i], ys[j], ts[k], nx, ny)] for i in 1:4, j in 1:4, k in 1:4]
        AFbyu = A_times_svec(F_times_svec(uvals))
        vvals = @SArray [Vs[_to_raw_index(xs[i], ys[j], ts[k], nx, ny)] for i in 1:4, j in 1:4, k in 1:4]
        AFbyv = A_times_svec(F_times_svec(vvals))

        T = eltype(AFbyu)
        res1 = zero(T)
        res2 = zero(T)
        for i in 1:4, j in 1:4, k in 1:4
            res1 += xp[i]*yp[j]*tp[k]*AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
            res2 += xp[i]*yp[j]*tp[k]*AFbyv[(i-1) + 4*(j-1) + 16*(k-1) + 1]
        end
    end
    return SVector{2}((res1/8.0, res2/8.0))
end

# TRICUBIC GRADIENT INTERPOLATION

function _get_cubic_grad_args(u, p, t)
    xs, x, ys, y, ts, t, nx, ny = _get_index2(u, p, t)

    xp = SVector{4}((1.0, x, x^2, x^3))
    dxp = SVector{4}((0.0, 1.0, 2*x, 3*x^2))
    yp = SVector{4}((1.0, y, y^2, y^3))
    dyp = SVector{4}((0.0, 1.0, 2*y, 3*y^2))
    tp = SVector{4}((1.0, t, t^2, t^3))

    return xs, ys, ts, xp, dxp, yp, dyp, tp, nx, ny
end

function _interp_grad_tricubic(Fs, xs, ys, ts, xp, dxp, yp, dyp, tp, nx, ny)
    T = promote_type(eltype(Fs), eltype(xp))
    result1 = zero(T)
    result2 = zero(T)
    result3 = zero(T)

    @inbounds begin
        fvals = @SArray [Fs[_to_raw_index(xs[i], ys[j], ts[k], nx, ny)] for i in 1:4, j in 1:4, k in 1:4]
        AFbyu = A_times_svec(F_times_svec(fvals))
        for i in 1:4, j in 1:4, k in 1:4
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
    end
    return SVector{3}((result1/8.0, result2/8.0, result3/8.0))
end

"""
    scalar_tricubic_gradient(u,p,t)

Calculates the (spatial) gradient of the interpolant obtained from [`scalar_tricubic`](@ref).
"""
function scalar_tricubic_gradient(u::StaticVector{2}, p::ItpMetadata, t)
    res = _interp_grad_tricubic(p.data[1], _get_cubic_grad_args(u, p, t)...)
    return SVector{2}((res[1]*p.nx/(p.UR[1] - p.LL[1]), res[2]*p.ny/(p.UR[2] - p.LL[2])))
end

"""
    uv_tricubic_eqvari(u,p,t)

The rhs for solving the linearized flow of the vector field (u,v) with `CoherentStructures.jl`.
"""
function uv_tricubic_eqvari(U::StaticMatrix{2,3}, p::ItpMetadata, t)
    u = SVector{2}((U[1,1], U[2,1]))
    args = _get_cubic_grad_args(u, p, t)    
    Uitp = _interp_grad_tricubic(p.data[1], args...)
    Vitp = _interp_grad_tricubic(p.data[2], args...)
    scale1 = p.nx / (p.UR[1] - p.LL[1])
    scale2 = p.ny / (p.UR[2] - p.LL[2])
    return @SMatrix [Uitp[3] (Uitp[1]*U[1,2]*scale1 + Uitp[2]*U[2,2]*scale2) (Uitp[1]*U[1,3]*scale1 + Uitp[2]*U[2,3]*scale2);
        Vitp[3] (Vitp[1]*U[1,2]*scale1 + Vitp[2]*U[2,2]*scale2) (Vitp[1]*U[1,3]*scale1 + Vitp[2]*U[2,3]*scale2)]
end

# interp_inbounds(x, y, p) = !isnan(p[3][x, y])

# inbounds_checker_bilinear(x, y, p) = !isnan(uv_trilinear(SVector{2,Float64}((x,y)), p[6], p[6][5][1])[1])

"""
    ssh_rhs(u, p, t)
Approximating geostrophic sea-surface velocities with the well-known formula

```math
u = -A(y)\\partial_y h(x,y,t)\\ v = A(y)\\partial_x h(x,y,t)
```

where:

* `u` -- longitudinal component of the velocity,
* `v` -- latitudinal component of the velocity,
* `x` -- longitude,
* `y` -- latitude,
* `h` -- sea-surface height.

and

```math
A(y) = g/(2 R^2 \\Omega \\sin y).
```

!!! note
    The sea surface height is assumed to be given in meters, the grid is assumed to
    be given in degrees (longitude/latitude). The output unit of the velocity is
    degree/day.
"""
function ssh_rhs(u, p::ItpMetadata, t)
    ∇h = scalar_tricubic_gradient(u, p, t)

    g = 9.807 #Gravitational constant (m/s^2)
    R = 6371e3 #Radius of earth (m)
    Ω = 7.2921159e-5*R #Mean angular velocity (m/s)

    # convert from m/deg to m/rad
    ∇h *= 360/(2π)

    #To convert ∇h to m/m, the first component must be divided by R*cos(deg2rad(u[2]))
    #and the second component must be divided by R, so we could do the R-division already.
    #But since later on we multiply by R again, we comment it out
    #∇h /= R 

    d2r = deg2rad(u[2])
    sind2r, cosd2r = sincos(d2r)

    #Calculate the velocity conversion factor in (ms)^{-1}
    C = g/(2*R*Ω*sind2r)

    #Now (except for the cos(deg2rad(u[2])) factor), C*skew(∇h) is in m/s,
    #but we want deg/s
    #We need to multiply the x-component by R*cos(deg2rad(u[2]) and the y-component by R
    #∇h *= R
    
    #We want deg/day instead of rad/s though
    C *= (24*3600)*360/(2π)

    return  SVector{2}((-∇h[2]*C/cosd2r, ∇h[1]*C*cosd2r))
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
"""
function uv_quadlinear(u::SVector{3,T}, p::ItpMetadata3{S}, t::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    Ws = p.data[3]
    return @inbounds _uv_quadlinear(u, Us, Vs, Ws, p, t)
end

function _uv_quadlinear(u::SVector{3,T}, Us::U, Vs::U, Ws::U, p::ItpMetadata3{S}, t::Float64) where {T<:Real,S,U}
    #Get data from p
    nx, ny, nz, nt = p.nx, p.ny, p.nz, p.nt
    ll1, ll2, ll3, t0 = p.LL
    ur1, ur2, ur3, tf = p.UR
    px, py, pt = p.periods

    xindex, xpp, xcoord = getIndex(u[1], ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ycoord = getIndex(u[2], ll2, ur2, py, ny, p.boundaryY)
    zindex, zpp, zcoord = getIndex(u[3], ll3, ur3, pz, nz, p.boundaryZ)
    tindex, tpp, tcoord = getIndex(t, t0, tf, pt, nt,  p.boundaryT)

    function singleTriLinear(what, zindex)
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
    return SVector{3,T}((singleQuadLinear(Us), singleQuadLinear(Vs), singleQuadLinear(Ws)))
end
