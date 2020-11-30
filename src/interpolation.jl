
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
            push!(coeff_updates.args, :($(Symbol("y", row)) += $val * $x[$col]))
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
number `x`, where `x` is in c_i, c_j and the interval [x0,xf) is partitioned into intervals
`[x0 = c_0, c_1), [c_1,c_2), ..., [c_(nx-1), xf)`, where all intervals have equal length.
"""
@inline function getIndex(x, x0, xf, xper, nx, boundary::BoundaryBehaviour)
    if boundary == periodic
        #spacing = (xf-x0)/nx
        #xindex, xcoord = gooddivrem(mod(x - x0, xf-x0),spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x - x0, (xf-x0))*(nx))/(xf-x0), 1)
        xpp = (xindex+1) % nx
    elseif boundary == flat
        #spacing = (xf-x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x-x0)*nx)/(xf-x0), 1)
        xpp = max(min(xindex + 1,nx-1),0)
        xindex = max(min(xindex,nx-1),0)
    elseif boundary  == outofbounds
        #spacing = (xf - x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x-x0)*nx)/(xf-x0), 1)
        xpp = xindex + 1
        if xpp >= nx || xindex < 0
                error("Out of bounds access at coordinate $x")
        end
    else #boundary == semiperiodic
        #spacing = mod(xf-x0,xper)/nx
        #xindex, xcoord = gooddivrem(mod(x-x0,xper), spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x-x0,xper)*nx)/mod(xf-x0,xper), 1)
        xpp = xindex + 1
        if xpp >= nx || xindex < 0
                error("Out of bounds access at coordinate $x")
        end
    end
    return xindex, xpp, xcoord
end

"""
    getIndex2(args)

Like `getIndex`, but also gives back one more set of points to the right/left.
"""
function getIndex2(x, x0, xf,xper, nx, boundary::BoundaryBehaviour)
    if boundary == periodic
        #spacing = (xf-x0)/nx
        #xindex, xcoord = gooddivrem(mod(x - x0, xf-x0),spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem((mod(x - x0, (xf - x0))*nx) / (xf - x0), 1)
        xindex = xindex % nx
        xpp = (xindex + 1) % nx
        xpp2 = (xindex + 2) % nx
        xmm = mod((xindex - 1), nx)
    elseif boundary == flat
        #spacing = (xf - x0)/nx
        #xindex, xcoord = gooddivrem(x-x0, spacing)
        #xcoord /= spacing
        xindex, xcoord = gooddivrem(((x - x0)*nx) / (xf - x0), 1)
        xindex = max(min(xindex, nx-1), 0)
        xpp = max(min(xindex + 1, nx-1), 0)
        xpp2 = max(min(xindex + 2, nx-1), 0)
        xmm = max(min(xindex - 1, nx-1), 0)
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
    return xindex, xpp, xpp2, xmm, xcoord
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
        data::T, boundaryX::BoundaryBehaviour, boundaryY::BoundaryBehaviour, bounradyZ::BoundaryBehaviour, boundaryT::BoundaryBehaviour) where {T}
    nx = length(xspan)
    ny = length(yspan)
    nz = length(zspan)
    nt = length(tspan)
    LL = SVector{4}((minimum(xspan), minimum(yspan), minimum(zspan), minimum(tspan)))
    UR = SVector{4}((maximum(xspan)+step(xspan), maximum(yspan)+step(yspan), maximum(zspan)+step(zspan), maximum(tspan)+step(tspan)))
    return ItpMetadata3(nx, ny, nz, nt, LL, UR, data, boundaryX, boundaryY, boundaryZ, boundaryT)
end

"""
    uv_trilinear(x, p, t)

Trilinear interpolation of velocity field at `x` at time `t`.
Velocity field stored in `p.data[1]` and `p.data[2]`.
"""
function uv_trilinear(x::SVector{2,<:Real}, p::ItpMetadata, t)
    #Get data from p
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px, py, pt = p.periods

    xindex, xpp, xcoord = getIndex(x[1], ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ycoord = getIndex(x[2], ll2, ur2, py, ny, p.boundaryY)
    tindex, tpp, tcoord = getIndex(t,    t0,  tf,  pt, nt, p.boundaryT)

    return SVector{2}(
        (
        _interp_trilinear(p.data[1], xindex, xpp, xcoord, yindex, ypp, ycoord, tindex, tpp, tcoord),
        _interp_trilinear(p.data[2], xindex, xpp, xcoord, yindex, ypp, ycoord, tindex, tpp, tcoord),
        )
    )
end

@inline function get_locations(xindex, yindex, tindex, xpp, ypp, tpp, xpp2, ypp2, tpp2, xmm, ymm, tmm, nx, ny, x::T, y::T, t::T) where {T}
    r = StaticArrays.SUnitRange(-1,2)
    locs = map(Iterators.product(r, r, r)) do (i, j, k)
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
        return xi + yi*nx + ti*nx*ny + 1
    end
    return locs
end

function _interp_tricubic(locs, xp, yp, tp, Fs)
    fvals = Fs[locs]
    AFbyu = A_times_svec(F_times_svec(fvals))

    res = zero(eltype(AFbyu))
    for i in 1:4, j in 1:4
        @simd for k in 1:4
            res += xp[i]*yp[j]*tp[k]*AFbyu[(i-1) + 4*(j-1) + 16*(k-1) + 1]
        end
    end

    return res/8.0
end

function base_tricubic_interpolation_gradient(
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

    function toRawIndex(i, j, k)
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
        return xi + yi*nx + ti*nx*ny + 1
    end

    @inbounds begin
        uvals = @SArray T[Us[toRawIndex(i,j,k)] for i in -1:2, j in -1:2, k in -1:2]
        AFbyu = A_times_svec(F_times_svec(uvals))
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
    return SVector{3,T}((result1/8.0, result2/8.0, result3/8.0))
end

"""
    uv_tricubic(x, p, t)

Component wise tricubic inuerpolation (Lekien-Marsden + finite differences for
values not specified in their paper) of velocity field at `x` at time `t`.
Velocity component fields are stored in `p.data[1]` and `p.data[2]`.
"""
function uv_tricubic(u::SVector{2,<:Real}, p::ItpMetadata, t)
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px, py, pt = p.periods

    xindex, xpp, xpp2, xmm, x = getIndex2(u[1], ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ypp2, ymm, y = getIndex2(u[2], ll2, ur2, py, ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, t = getIndex2(t,    t0,  tf,  pt, nt, p.boundaryT)

    locs = get_locations(xindex, yindex, tindex, xpp, ypp, tpp, xpp2, ypp2, tpp2, xmm, ymm, tmm, nx, ny, x, y, t)
    xp = SVector{4}((1.0, x, x^2, x^3))
    yp = SVector{4}((1.0, y, y^2, y^3))
    tp = SVector{4}((1.0, t, t^2, t^3))

    return SVector{2}(_interp_tricubic(locs, xp, yp, tp, p.data[1]), _interp_tricubic(locs, xp, yp, tp, p.data[2]))
    
    # base_tricubic_interpolation(
    #     xindex, yindex, tindex,
    #     xpp, ypp, tpp,
    #     xpp2, ypp2, tpp2,
    #     xmm, ymm, tmm,
    #     nx, ny,
    #     xcoord, ycoord, tcoord,
    #     p.data[1], p.data[2],
    # )
end

"""
    scalar_trilinear(x, p, t)

Trilinear interpolation of scalar field at `x` at time `t`.
Scalar field is assumed to be stored in `p.data[1]`.
"""
function scalar_trilinear(x::SVector{2}, p, t::Float64)
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px, py, pt = p.periods

    xindex, xpp, xcoord = getIndex(x[1], ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ycoord = getIndex(x[2], ll2, ur2, py, ny, p.boundaryY)
    tindex, tpp, tcoord = getIndex(t,    t0,  tf,  pt, nt, p.boundaryT)

    return @inbounds _interp_trilinear(p.data[1], xindex, xpp, xcoord, yindex, ypp, ycoord, tindex, tpp, tcoord)
end

@inline function _interp_trilinear(Fs, xindex, xpp, xcoord, yindex, ypp, ycoord, tindex, tpp, tcoord)
    r1u = Fs[xindex + 1, yindex + 1, tindex + 1]*(1-xcoord) + Fs[xpp + 1, yindex + 1, tindex + 1]*xcoord
    r2u = Fs[xindex + 1, ypp + 1, tindex + 1]*(1-xcoord) + Fs[xpp + 1, ypp + 1, tindex + 1]*xcoord
    r3u = Fs[xindex + 1, yindex + 1, tpp + 1]*(1-xcoord) + Fs[xpp + 1, yindex + 1, tpp + 1]*xcoord
    r4u = Fs[xindex + 1, ypp + 1, tpp + 1]*(1-xcoord) + Fs[xpp + 1,ypp + 1, tpp + 1 ]*xcoord
    return ((1-tcoord)*((1-ycoord)*r1u + ycoord*r2u) + tcoord*((1-ycoord)*r3u + ycoord*r4u))
end

"""
    scalar_tricubic(x, p, t)

Tricubic interpolation (Lekien-Marsden + finite differences for values not
specified in their paper) of scalar field at `x` at time `t`.
Scalar field is assumed to be stored in `p.data[1]`.
"""
function scalar_tricubic(u::StaticVector{2,<:Real}, p, t::Float64)
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px, py, pt = p.periods

    xindex, xpp, xpp2, xmm, x = getIndex2(u[1], ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ypp2, ymm, y = getIndex2(u[2], ll2, ur2, py, ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, t = getIndex2(t,    t0,  tf,  pt, nt, p.boundaryT)

    locs = get_locations(xindex, yindex, tindex, xpp, ypp, tpp, xpp2, ypp2, tpp2, xmm, ymm, tmm, nx, ny, x, y, t)
    xp = SVector{4}((1.0, x, x^2, x^3))
    yp = SVector{4}((1.0, y, y^2, y^3))
    tp = SVector{4}((1.0, t, t^2, t^3))
    
    return _interp_tricubic(locs, xp, yp, tp, p.data[1])
end

"""
    scalar_tricubic_gradient(u,p,t)

Calculates the (spatial) gradient of the function used in scalar_tricubic
"""
function scalar_tricubic_gradient(u::StaticVector{2,T}, p::ItpMetadata{S}, t::Float64) where {T<:Real,S}
    sshs = p.data[1]
    return @inbounds _scalar_tricubic_gradient(u, sshs, p, t)
end


function _scalar_tricubic_gradient(u::StaticVector{2,T}, sshs::U, p::ItpMetadata{S}, t::Float64) where {T<:Real,S,U}
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR

    xindex, xpp, xpp2, xmm, xcoord = getIndex2(u[1],ll1,ur1, p.periods[1], nx, p.boundaryX)
    yindex, ypp, ypp2, ymm, ycoord = getIndex2(u[2],ll2,ur2, p.periods[2], ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, tcoord = getIndex2(t, t0, tf, p.periods[3], nt, p.boundaryT)

    res1 = base_tricubic_interpolation_gradient(xindex, yindex, tindex,
                                                xpp, ypp, tpp,
                                                xpp2, ypp2, tpp2,
                                                xmm, ymm, tmm,
                                                nx, ny,
                                                xcoord, ycoord, tcoord,
                                                sshs)

    return SVector{2,T}((res1[1]*nx/(ur1-ll1), res1[2]*ny/(ur2-ll2)))
end


"""
    uv_tricubic_eqvari(u,p,t)

The rhs for solving the linearized flow of the vector field (u,v) with CoherentStructures.jl
"""
function uv_tricubic_eqvari(u::StaticMatrix{2,3,T}, p::ItpMetadata{S}, t::Float64) where {T<:Real,S}
    Us = p.data[1]
    Vs = p.data[2]
    return @inbounds _uv_tricubic_eqvari(u, Us, Vs, p, t)
end

function _uv_tricubic_eqvari(uIn::StaticMatrix{2,3}, Us::U, Vs::U, p::ItpMetadata, t::Float64) where {U}
    u = SVector{2}((uIn[1,1], uIn[2,1]))
    nx, ny, nt = p.nx, p.ny, p.nt
    ll1, ll2, t0 = p.LL
    ur1, ur2, tf = p.UR
    px,py,pt = p.periods

    xindex, xpp, xpp2, xmm, xcoord = getIndex2(u[1], ll1, ur1, px, nx, p.boundaryX)
    yindex, ypp, ypp2, ymm, ycoord = getIndex2(u[2], ll2, ur2, py, ny, p.boundaryY)
    tindex, tpp, tpp2, tmm, tcoord = getIndex2(t,    t0,  tf,  pt, nt, p.boundaryT)

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
    Vitp[3] (Vitp[1]*uIn[1,2]*nx/(ur1-ll1) + Vitp[2]*uIn[2,2]*ny/(ur2-ll2)) (Vitp[1]*uIn[1,3]*nx/(ur1-ll1) + Vitp[2]*uIn[2,3]*ny/(ur2-ll2))]
end

function interp_inbounds(x, y, p)
    return !isnan(p[3][x, y])
end
function inbounds_checker_bilinear(x, y, p)
    return !isnan(uv_trilinear(SVector{2,Float64}((x,y)), p[6], p[6][5][1])[1])
end

"""
    ssh_rhs(u, p, t)

Approximating geostrophic sea-surface velocities with the well-known formula

```math
u = -A(y)\\partial_y h(x,y,t),\\ v = A(y)\\partial_x h(x,y,t),
```

where:

* `u` -- longitudinal component of the velocity,
* `v` -- latitudinal component of the velocity,
* `x` -- longitude,
* `y` -- latitude,
* `h` -- sea-surface height,

and

```math
A(y) = g/(R^2 2 \\Omega \\sin y).
```
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

    return  SVector{2,Float64}((-∇h[2]*C/cosd2r, ∇h[1]*C*cosd2r))
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
