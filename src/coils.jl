abstract type AbstractCoil{T<:Real,C<:Real} end

mutable struct PointCoil{T<:Real,C<:Real} <: AbstractCoil{T,C}
    R::T
    Z::T
    current::C
end

PointCoil(R, Z) = PointCoil(R, Z, 0.0)

"""
    ParallelogramCoil{T,C} <:  AbstractCoil{T,C}

Parallelogram coil with the R, Z, ΔR, ΔZ, θ₁, θ₂ formalism (as used by EFIT, for example)
Here θ₁ and θ₂ are the shear angles along the x- and y-axes, respectively, in degrees.
"""
mutable struct ParallelogramCoil{T<:Real,C<:Real} <: AbstractCoil{T,C}
    R::T
    Z::T
    ΔR::T
    ΔZ::T
    θ₁::T
    θ₂::T
    spacing::T
    current::C
end

ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing) = ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing, 0.0)

function ParallelogramCoil(R::T, Z::T, ΔR::T, ΔZ::T, θ₁::T, θ₂::T; spacing::T=0.01) where {T<:Real}
    return ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing)
end

mutable struct DistributedCoil{T<:Real,C<:Real} <: AbstractCoil{T,C}
    R::Vector{T}
    Z::Vector{T}
    current::C
end

DistributedCoil(R, Z) = DistributedCoil(R, Z, 0.0)
"""
    DistributedParallelogramCoil(Rc::T, Zc::T, ΔR::T, ΔZ::T, θ₁::T, θ₂::T; spacing::T=0.01) where {T<:Real}

NOTE: if spacing <= 0.0 then current filaments are placed at the vertices
"""
function DistributedParallelogramCoil(Rc::T, Zc::T, ΔR::T, ΔZ::T, θ₁::T, θ₂::T; spacing::T=0.01) where {T<:Real}
    if spacing <= 0.0
        dR = [-0.5 * ΔR, 0.5 * ΔR]
        dZ = [-0.5 * ΔZ, 0.5 * ΔZ]
    else
        dR = LinRange(-0.5 * ΔR, 0.5 * ΔR, Int(ceil(1.0 + ΔR / spacing)))
        dZ = LinRange(-0.5 * ΔZ, 0.5 * ΔZ, Int(ceil(1.0 + ΔZ / spacing)))
    end

    α₁ = tan(π * θ₁ / 180.0)
    α₂ = tan(π * (θ₂ + 90.0) / 180.0)

    Nr = length(dR)
    Nz = length(dZ)
    N = Nr * Nz
    Z = Array{Float64,1}(undef, N)
    R = Array{Float64,1}(undef, N)
    for i in 1:Nr
        for j in 1:Nz
            k = j + (i - 1) * Nz
            @inbounds R[k] = Rc + dR[i] - α₂ * dZ[j]
            @inbounds Z[k] = Zc + dZ[j] + α₁ * dR[i]
        end
    end
    return DistributedCoil(R, Z)
end

function DistributedCoil(C::ParallelogramCoil)
    return DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ₁, C.θ₂; C.spacing)
end

# ================ #
# encircling_coils #
# ================ #

function encircling_coils(EQfixed::MXHEquilibrium.AbstractEquilibrium, n_coils::Integer)
    bnd = MXHEquilibrium.plasma_boundary(EQfixed; precision=0.0)
    if bnd === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        bnd = MXHEquilibrium.plasma_boundary(EQfixed)
    end

    return encircling_coils(bnd.r, bnd.z, n_coils)
end

function encircling_coils(shot::TEQUILA.Shot, n_coils::Integer)
    bnd_r, bnd_z = MillerExtendedHarmonic.MXH(shot.surfaces[:, end])()
    return encircling_coils(bnd_r, bnd_z, n_coils)
end

function encircling_coils(bnd_r::AbstractVector{T}, bnd_z::AbstractVector{T}, n_coils::Integer) where {T<:Real}
    mxh = MillerExtendedHarmonic.MXH(bnd_r, bnd_z, 2)
    mxh.R0 = mxh.R0 .* (1.0 + mxh.ϵ)
    mxh.ϵ = 0.9
    Θ = LinRange(0, 2π, n_coils + 1)[1:end-1]
    return [PointCoil(r, z) for (r, z) in mxh.(Θ)]
end

# ============== #
#   Convex hull  #
# ============== #

@inline function points_isless(p::AbstractVector{T}, q::AbstractVector{T}) where {T}
    return p[1] < q[1] || (p[1] == q[1] && p[2] < q[2])
end

@inline function points_isless(p::Tuple{T,T}, q::Tuple{T,T}) where {T}
    return p[1] < q[1] || (p[1] == q[1] && p[2] < q[2])
end

@inline function isrightturn(p::AbstractVector{T}, q::AbstractVector{T}, r::AbstractVector{T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) < 0.0
end

@inline function isrightturn(p::Tuple{T,T}, q::Tuple{T,T}, r::Tuple{T,T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) < 0.0
end

function halfhull(points::AbstractVector)
    halfhull = similar(points)
    n = 0
    for p in points
        while n > 1 && !isrightturn(halfhull[n-1], halfhull[n], p)
            n -= 1
        end
        n += 1
        halfhull[n] = p
    end
    return view(halfhull, 1:n)
end

function grahamscan!(points::AbstractVector)
    sort!(points; lt=points_isless)
    upperhull = halfhull(points)
    reverse!(points)
    lowerhull = halfhull(points)
    return [upperhull; lowerhull[2:end-1]]
end

function convex_hull!(xy_points::AbstractVector; closed_polygon::Bool)
    hull = grahamscan!(xy_points)
    if closed_polygon
        return push!(hull, hull[1])
    else
        return hull
    end
end

function convex_hull(xy_points::AbstractVector; closed_polygon::Bool)
    return convex_hull!(deepcopy(xy_points); closed_polygon)
end

function convex_hull(x::AbstractVector{T}, y::AbstractVector{T}; closed_polygon::Bool) where {T}
    xy_points = [(xx, yx) for (xx, yx) in zip(x, y)]
    return convex_hull!(xy_points; closed_polygon)
end

function convex_hull(C::ParallelogramCoil)
    C = DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ₁, C.θ₂; spacing=0.0)
    return collect(zip(C.R, C.Z))
end

function convex_hull(C::DistributedCoil)
    return convex_hull(C.R, C.Z; closed_polygon=true)
end

# ======== #
#   Plot   #
# ======== #
@recipe function plot_coil(C::ParallelogramCoil)
    DistributedCoil(C)
end

@recipe function plot_coil(C::DistributedCoil)
    hull = convex_hull(C)
    R = [r for (r, z) in hull]
    Z = [z for (r, z) in hull]
    @series begin
        color --> :black
        alpha --> 0.2
        label --> ""
        Shape(R,Z)
    end
end

@recipe function plot_coil(C::PointCoil)
    @series begin
        seriestype --> :scatter
        marker --> :circle
        markercolor --> :black
        [C.R], [C.Z]
    end
end