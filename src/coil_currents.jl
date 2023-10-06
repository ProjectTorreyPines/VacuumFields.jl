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

# ============== #
#   Convex hull  #
# ============== #
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
function plot_coil(C::ParallelogramCoil)
    return plot_coil(DistributedCoil(C))
end

function plot_coil(C::DistributedCoil)
    hull = convex_hull(C)
    R = [r for (r, z) in hull]
    Z = [z for (r, z) in hull]
    return plot!(Shape(R, Z); color=:black, alpha=0.2, label="")
end

function plot_coil(C::PointCoil)
    return plot!([C.R], [C.Z]; marker=:circle, markercolor=:black)
end

# ======== #
#   Green  #
# ======== #
function Green(C::ParallelogramCoil, R::T, Z::T, scale_factor::Float64=1.0) where {T<:Real}
    return Green(DistributedCoil(C), R, Z, scale_factor)
end

function Green(C::DistributedCoil, R::T, Z::T, scale_factor::Float64=1.0) where {T<:Real}
    return sum(Green(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end

function Green(C::PointCoil, R::T, Z::T, scale_factor::Float64=1.0) where {T<:Real}
    return Green(C.R, C.Z, R, Z, scale_factor)
end

@inline function Green(X::T, Y::T, R::T, Z::T, scale_factor::Float64=1.0) where {T<:Real}
    XR = X * R
    m = 4.0 * XR / ((X + R)^2 + (Y - Z)^2) # this is k^2
    Km, Em = ellipke(m)
    return inv2π * (2.0 * Em - (2.0 - m) * Km) * sqrt(XR / m) * scale_factor
end

# ======== #
#   Utils  #
# ======== #
function cumlength(R::T, Z::T) where {T<:AbstractVector{Float64}}
    # Length along boundary
    L = Array{Float64,1}(undef, length(R))
    @inbounds L[1] = 0
    for i in eachindex(R)[2:end]
        @inbounds L[i] = L[i-1] + sqrt((R[i] - R[i-1])^2 + (Z[i] - Z[i-1])^2)
    end
    return L
end

function bounds(Cs::AbstractVector{<:AbstractCoil})
    Rmin = minimum(Cs[1].R)
    Rmax = maximum(Cs[1].R)
    Zmin = minimum(Cs[1].Z)
    Zmax = maximum(Cs[1].Z)
    for C in Cs
        Rmin = min(Rmin, minimum(C.R))
        Rmax = max(Rmax, maximum(C.R))
        Zmin = min(Zmin, minimum(C.Z))
        Zmax = max(Zmax, maximum(C.Z))
    end
    return Rmin, Rmax, Zmin, Zmax
end

# ========== #
#   Physics  #
# ========== #
function fixed_boundary(EQfixed::MXHEquilibrium.AbstractEquilibrium, Sb::MXHEquilibrium.PlasmaBoundary)
    Rb, Zb = Sb.r, Sb.z
    Lb = cumlength(Rb, Zb)

    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    fac = EQfixed.cocos.sigma_rhotp * EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    dψdn_R = fac .* MXHEquilibrium.poloidal_Bfield.(Ref(EQfixed), Rb, Zb)
    return (Rb, Zb, Lb, dψdn_R)
end

function fixed_boundary(shot::TEQUILA.Shot, bnd::MillerExtendedHarmonic.MXH; Nb::Integer=100 * length(bnd.c))
    return fixed_boundary(shot, bnd(Nb; adaptive=false)...)
end

function fixed_boundary(shot::TEQUILA.Shot, Rb::AbstractVector{T}, Zb::AbstractVector{T}) where {T<:Real}
    Lb = cumlength(Rb, Zb)
    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    cocos = MXHEquilibrium.cocos(shot)
    fac = cocos.sigma_rhotp * cocos.sigma_Bp * (2π)^cocos.exp_Bp
    dψdn_R = fac .* MXHEquilibrium.poloidal_Bfield.(Ref(shot), Rb, Zb)
    return (Rb, Zb, Lb, dψdn_R)
end

"""
    ψp_on_fixed_eq_boundary(EQfixed,
                            fixed_coils=[],
                            ψbound=0.0;
                            Rx=[], Zx=[],
                            fraction_inside=0.999)

Calculate ψ from image currents on boundary at surface p near boundary
"""
function ψp_on_fixed_eq_boundary(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    fixed_coils::Vector{<:AbstractCoil{T,C}}=PointCoil{T,C}[],
    ψbound::Real=0.0;
    Rx::AbstractVector{Float64}=Float64[],
    Zx::AbstractVector{Float64}=Float64[],
    fraction_inside::Float64) where {T<:Real,C<:Real}

    ψ0, ψb = MXHEquilibrium.psi_limits(EQfixed)
    ψb, Sb = MXHEquilibrium.plasma_boundary_psi(EQfixed; precision=0.0)
    if Sb === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        ψb, Sb = MXHEquilibrium.plasma_boundary_psi(EQfixed)
    end

    # ψp is the flux from the image currents on the plasma boundary
    # which is equal and opposite to the flux from the plasma current
    # ψ₀ = ψp + ψim
    # ψ₁ = ψp + ψcoil
    # ψcoil = (ψ₁ - ψ₀) + ψim
    Sp = MXHEquilibrium.flux_surface(EQfixed, fraction_inside * (ψb - ψ0) + ψ0)
    n = Int(floor(length(Sp.r) / 100.0)) + 1 # roughly at most 100 points
    Rp, Zp = Sp.r[1:n:end-1], Sp.z[1:n:end-1]

    append!(Rp, Rx)
    append!(Zp, Zx)

    # this is the image current contribution to the control points
    ψp = Array{Float64,1}(undef, length(Rp))
    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed, Sb)
    @threads for i in eachindex(Rp)
        ψp[i] = -trapz(Lb, dψdn_R .* Green.(Rb, Zb, Rp[i], Zp[i]))
    end

    # this is to account that the control points are inside of the LCFS
    # applied only to the plasma points
    ψp[1:end-length(Rx)] .-= (1.0 - fraction_inside) * (ψb - ψ0)

    # add in desired boundary flux
    ψbound != 0.0 && (ψp .+= ψbound)

    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp

    # Compute flux from fixed coils and subtract from ψp to match
    # This works whether ψp is a constant or a vector
    ψfixed = zeros(length(Rp))
    @threads for j in eachindex(fixed_coils)
        for i in eachindex(Rp)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rp[i], Zp[i]) * fixed_coils[j].current
        end
    end
    ψp = ψp .- ψfixed

    return Bp_fac, ψp, Rp, Zp
end

# Compute the plasma contribution to ψ at points fraction_inside the shot boundary,
#   optionally including vacuum points given by Rx and Zx
function ψp_on_fixed_eq_boundary(
    shot::TEQUILA.Shot,
    fixed_coils::Vector{<:AbstractCoil{T,C}}=PointCoil{T,C}[],
    ψbound::Real=0.0;
    Rx::AbstractVector{Float64}=Float64[],
    Zx::AbstractVector{Float64}=Float64[],
    fraction_inside::Float64=0.999) where {T<:Real,C<:Real}

    Sb = @views MillerExtendedHarmonic.MXH(shot.surfaces[:, end])

    # ψp is the flux from the image currents on the plasma boundary
    # which is equal and opposite to the flux from the plasma current
    # ψ₀ = ψp + ψim
    # ψ₁ = ψp + ψcoil
    # ψcoil = (ψ₁ - ψ₀) + ψim
    Sp = MillerExtendedHarmonic.MXH(shot, fraction_inside)
    Np = 99
    Nx = length(Rx)

    Rp = zeros(Np + Nx)
    Zp = zeros(Np + Nx)
    r, z = Sp(Np + 1; adaptive=false)
    Rp[1:Np] .= r[1:end-1]
    Zp[1:Np] .= z[1:end-1]
    if Nx > 0
        Rp[(Np+1):end] .= Rx
        Zp[(Np+1):end] .= Zx
    end

    Rb, Zb, Lb, dψdn_R = @views fixed_boundary(shot, Sb)

    # this is the image current contribution to the control points
    ψp = Array{Float64,1}(undef, length(Rp))
    @threads for i in eachindex(Rp)
        ψp[i] = -trapz(Lb, dψdn_R .* Green.(Rb, Zb, Rp[i], Zp[i]))
    end

    # this is to account that the control points are inside of the LCFS
    # applied only to the plasma points
    ψp[1:Np] .+= TEQUILA.psi_ρθ(shot, fraction_inside, 0.0)

    # add in desired boundary flux
    ψbound != 0.0 && (ψp .+= ψbound)

    cocos = MXHEquilibrium.cocos(shot)
    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp

    # Compute flux from fixed coils and subtract from ψp to match
    # This works whether ψp is a constant or a vector
    ψfixed = zeros(length(Rp))
    @threads for j in eachindex(fixed_coils)
        for i in eachindex(Rp)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rp[i], Zp[i]) * fixed_coils[j].current
        end
    end
    ψp = ψp .- ψfixed

    return Bp_fac, ψp, Rp, Zp
end

# Compute the plasma contribution to ψ at control points given by Rc and Zc
function ψp_on_fixed_eq_boundary(
    shot::TEQUILA.Shot,
    Rc::AbstractVector{Float64},
    Zc::AbstractVector{Float64};
    fixed_coils::Vector{<:AbstractCoil{T,C}}=PointCoil{T,C}[],
    ψbound::Real=0.0) where {T<:Real,C<:Real}

    Sb = @views MillerExtendedHarmonic.MXH(shot.surfaces[:, end])

    Rb, Zb, Lb, dψdn_R = @views fixed_boundary(shot, Sb)

    # this is the image current contribution to the control points
    ψp = Array{Float64,1}(undef, length(Rc))
    @threads for i in eachindex(Rc)
        ψp[i] = -trapz(Lb, dψdn_R .* Green.(Rb, Zb, Rc[i], Zc[i]))
    end

    ψp .+= shot.(Rc, Zc)

    # add in desired boundary flux
    ψbound != 0.0 && (ψp .+= ψbound)

    cocos = MXHEquilibrium.cocos(shot)
    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp

    # Compute flux from fixed coils and subtract from ψp to match
    # This works whether ψp is a constant or a vector
    ψfixed = zeros(length(Rc))
    @threads for j in eachindex(fixed_coils)
        for i in eachindex(Rc)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rc[i], Zc[i]) * fixed_coils[j].current
        end
    end
    ψp = ψp .- ψfixed

    return Bp_fac, ψp, Rc, Zc
end

"""
    field_null_on_boundary(
        ψp_constant::Real,
        Rp::AbstractVector{Float64},
        Zp::AbstractVector{Float64},
        fixed_coils::Vector{<:AbstractCoil{T,C}}=PointCoil{T,C}[],
        ψbound::Real=0.0,
        cocos::Int=11) where {T<:Real,C<:Real}

Account for effect of fixed coils on ψp_constant
"""
function field_null_on_boundary(
    ψp_constant::Real,
    Rp::AbstractVector{Float64},
    Zp::AbstractVector{Float64},
    fixed_coils::Vector{<:AbstractCoil{T,C}}=PointCoil{T,C}[],
    ψbound::Real=0.0,
    cocos::Int=11) where {T<:Real,C<:Real}

    # add in desired boundary flux
    ψbound != 0.0 && (ψp_constant .+= ψbound)

    Bp_fac = MXHEquilibrium.cocos(cocos).sigma_Bp * (2π)^MXHEquilibrium.cocos(cocos).exp_Bp

    # Compute flux from fixed coils and subtract from ψp to match
    ψfixed = zeros(length(Rp))
    @threads for j in eachindex(fixed_coils)
        for i in eachindex(Rp)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rp[i], Zp[i]) * fixed_coils[j].current
        end
    end
    ψp = ψp_constant .- ψfixed

    return Bp_fac, ψp, Rp, Zp
end

function currents_to_match_ψp(
    Bp_fac::Real,
    ψp::AbstractVector{<:Real},
    Rp::AbstractVector{Float64},
    Zp::AbstractVector{Float64},
    coils::AbstractVector{<:AbstractCoil{T,C}};
    weights::Vector{Float64}=Float64[],
    λ_regularize::Float64=1E-16,
    return_cost::Bool=false) where {T<:Real,C<:Real}

    # Compute coil currents needed to recreate ψp at points (Rp,Zp)
    # Build matrix relating coil Green's functions to boundary points
    Gcp = Array{T,2}(undef, length(Rp), length(coils))
    Threads.@threads for j in eachindex(coils)
        for i in eachindex(Rp)
            @inbounds Gcp[i, j] = μ₀ * Bp_fac * Green(coils[j], Rp[i], Zp[i])
        end
    end

    # handle weights
    if length(weights) > 0
        ψp = weights .* ψp
        mul!(Gcp, Diagonal(weights), Gcp)
    end

    # Least-squares solve for coil currents
    if λ_regularize > 0
        # Least-squares with regularization
        # https://www.youtube.com/watch?v=9BckeGN0sF0
        reg_solve(A, b, λ) = inv(A' * A + λ * I) * A' * b
        Ic0 = reg_solve(Gcp, ψp, λ_regularize / length(coils)^2)
    else
        Ic0 = Gcp \ ψp
    end

    # update values of coils current
    for k in eachindex(coils)
        coils[k].current = Ic0[k]
    end

    if return_cost
        cost(Ic) = norm(Gcp * Ic .- ψp) / norm(ψp)
        return Ic0, cost(Ic0)
    else
        return Ic0
    end
end

function fixed_eq_currents(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:AbstractCoil{T,C}},
    fixed_coils::Vector{<:AbstractCoil{T,C}}=PointCoil{T,C}[],
    ψbound::Real=0.0;
    λ_regularize::Float64=1E-16,
    return_cost::Bool=false,
    fraction_inside::Float64=1.0 - 1E-6) where {T<:Real,C<:Real}

    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(EQfixed, fixed_coils, ψbound)
    return currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; fraction_inside, λ_regularize, return_cost)
end


# ***************************************************
# Transform fixed-boundary ψ to free-boundary ψ
# ***************************************************
"""
    fixed2free(
        EQfixed::MXHEquilibrium.AbstractEquilibrium,
        n_coils::Integer;
        Rx::AbstractVector{Float64}=Float64[],
        Zx::AbstractVector{Float64}=Float64[],
        fraction_inside::Float64=1.0 - 1E-6,
        λ_regularize::Float64=0.0,
        Rgrid::AbstractVector{Float64}=EQfixed.r,
        Zgrid::AbstractVector{Float64}=EQfixed.z)

Distribute n point coils around fixed boundary plasma to get a free boundary ψ map
"""
function fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    n_coils::Integer;
    Rx::AbstractVector{Float64}=Float64[],
    Zx::AbstractVector{Float64}=Float64[],
    fraction_inside::Float64=1.0 - 1E-6,
    λ_regularize::Float64=0.0,
    Rgrid::AbstractVector{Float64}=EQfixed.r,
    Zgrid::AbstractVector{Float64}=EQfixed.z)

    ψbound = MXHEquilibrium.psi_boundary(EQfixed; precision=0.0)
    if ψbound === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        ψbound = MXHEquilibrium.psi_boundary(EQfixed)
    end

    coils = encircling_coils(EQfixed, n_coils)
    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(EQfixed, coils, ψbound; fraction_inside, Rx, Zx)

    λ_regularize = optimal_λ_regularize(λ_regularize, Bp_fac, ψp, Rp, Zp, coils)
    currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; λ_regularize)
    return transpose(fixed2free(EQfixed, coils, Rgrid, Zgrid))
end

function cost_λ_regularize(
    λ_reg::Float64,
    Bp_fac::Real,
    ψp::AbstractVector{<:Real},
    Rp::AbstractVector{Float64},
    Zp::AbstractVector{Float64},
    coils::AbstractVector{<:AbstractCoil{T,C}}
) where {T<:Real,C<:Real}
    c = currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; λ_regularize=10^(λ_reg), return_cost=true)[2]
    return c^2
end

function optimal_λ_regularize(
    λ_regularize::Float64,
    Bp_fac::Real,
    ψp::AbstractVector{<:Real},
    Rp::AbstractVector{Float64},
    Zp::AbstractVector{Float64},
    coils::AbstractVector{<:AbstractCoil{T,C}}
) where {T<:Real,C<:Real}
    if λ_regularize == 0.0
        λ_range_exp = collect(-20:0.5:-10)
        cost_λ = [cost_λ_regularize(λ, Bp_fac, ψp, Rp, Zp, coils) for λ in λ_range_exp]
        λ_regularize = 10^λ_range_exp[argmin(cost_λ)]
    end
    return λ_regularize
end

function fixed2free(shot::TEQUILA.Shot, n_coils::Integer; n_grid=129, scale::Float64=1.2, kwargs...)
    R0 = shot.surfaces[1, end]
    Z0 = shot.surfaces[2, end]
    ϵ = shot.surfaces[3, end]
    κ = shot.surfaces[4, end]
    a = R0 * ϵ * scale
    b = κ * a

    Rgrid = range(max(R0 - a, 0.0), R0 + a, n_grid)
    Zgrid = range(Z0 - b, Z0 + b, n_grid)
    return Rgrid, Zgrid, fixed2free(shot, n_coils, Rgrid, Zgrid; kwargs...)
end

# On (Rgrid, Zgrid), convert the fixed equilibrium ψ to free using n_coils encircling coils
#   whose currents are determined by matching the ψ on shot's boundary and at (Rx, Zx)
function fixed2free(
    shot::TEQUILA.Shot,
    n_coils::Integer,
    Rgrid::AbstractVector{Float64},
    Zgrid::AbstractVector{Float64};
    Rx::AbstractVector{Float64}=Float64[],
    Zx::AbstractVector{Float64}=Float64[],
    fraction_inside::Union{Nothing,Float64}=1.0 - 1E-6,
    λ_regularize::Float64=0.0,
    ψbound::Real=0.0)

    coils = encircling_coils(shot, n_coils)
    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(shot, coils, ψbound; fraction_inside, Rx, Zx)

    λ_regularize = optimal_λ_regularize(λ_regularize, Bp_fac, ψp, Rp, Zp, coils)
    currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; λ_regularize)
    return transpose(fixed2free(shot, coils, Rgrid, Zgrid; ψbound))
end

# Convert the fixed equilibrium ψ to free using n_coils encircling coils on a bounding box around shot
#   whose currents are determined by matching the ψ at (Rc, Zc)
function fixed2free(
    shot::TEQUILA.Shot,
    n_coils::Integer,
    Rc::AbstractVector{<:Real},
    Zc::AbstractVector{<:Real};
    n_grid=129,
    scale::Float64=1.2,
    kwargs...)

    R0 = shot.surfaces[1, end]
    Z0 = shot.surfaces[2, end]
    ϵ = shot.surfaces[3, end]
    κ = shot.surfaces[4, end]
    a = R0 * ϵ * scale
    b = κ * a

    Rgrid = range(max(R0 - a, 0.0), R0 + a, n_grid)
    Zgrid = range(Z0 - b, Z0 + b, n_grid)
    return Rgrid, Zgrid, fixed2free(shot, n_coils, Rc, Zc, Rgrid, Zgrid; kwargs...)
end

# On (Rgrid, Zgrid), convert the fixed equilibrium ψ to free using n_coils encircling coils,
#   whose currents are determined by matching the ψ at (Rc, Zc)
function fixed2free(
    shot::TEQUILA.Shot,
    n_coils::Integer,
    Rc::AbstractVector{Float64},
    Zc::AbstractVector{Float64},
    Rgrid::AbstractVector{Float64},
    Zgrid::AbstractVector{Float64};
    λ_regularize::Float64=0.0,
    ψbound::Real=0.0)

    coils = encircling_coils(shot, n_coils)
    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(shot, Rc, Zc; ψbound)

    λ_regularize = optimal_λ_regularize(λ_regularize, Bp_fac, ψp, Rp, Zp, coils)
    weights = ones(length(Rp))
    weights[(end-11):end] .= 100.
    currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; λ_regularize)
    return transpose(fixed2free(shot, coils, Rgrid, Zgrid; ψbound))
end

function fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:ParallelogramCoil},
    R::AbstractVector{T},
    Z::AbstractVector{T}) where {T<:Real}
    dcoils = DistributedCoil.(coils)
    return fixed2free(EQfixed, dcoils, R, Z)
end

function fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:AbstractCoil{T,C}},
    R::AbstractVector{T},
    Z::AbstractVector{T}) where {T<:Real,C<:Real}

    # upsample tracing of boundary to get smooth LCFS from fixed2free
    R1 = LinRange(minimum(R), maximum(R), length(R) * 10)
    Z1 = LinRange(minimum(Z), maximum(Z), length(Z) * 10)
    ψb, Sb = MXHEquilibrium.plasma_boundary_psi(EQfixed; precision=0.0, r=R1, z=Z1)
    if Sb === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        ψb, Sb = MXHEquilibrium.plasma_boundary_psi(EQfixed; r=R1, z=Z1)
    end

    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    ψ_f2f = T[MXHEquilibrium.in_boundary(Sb, (r, z)) ? EQfixed(r, z) : ψb for z in Z, r in R]

    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed, Sb)

    # ψ from image and coil currents
    Vbs = [similar(Lb) for _ in 1:Threads.nthreads()] # preallocate
    Threads.@threads for (k, indexes) in collect(enumerate(chunk_indices((length(R), length(Z)), Threads.nthreads())))
        @inbounds for index in indexes
            i = index[1]
            j = index[2]
            r = R[i]
            z = Z[j]
            Vb = Vbs[k]
            Vb .= dψdn_R .* Green.(Rb, Zb, r, z)
            ψi = -Trapz.trapz(Lb, Vb)
            ψc = μ₀ * Bp_fac * sum(coil.current * Green(coil, r, z) for coil in coils)
            ψ_f2f[j, i] += ψc - ψi
        end
    end

    return ψ_f2f
end

function fixed2free(
    shot::TEQUILA.Shot,
    coils::AbstractVector{<:AbstractCoil{T,C}},
    R::AbstractVector{T},
    Z::AbstractVector{T};
    ψbound::Real=0.0) where {T<:Real,C<:Real}

    bnd = @views MillerExtendedHarmonic.MXH(shot.surfaces[:, end])
    cocos = MXHEquilibrium.cocos(shot)
    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp
    ψ_f2f = T[shot(r, z) + ψbound for z in Z, r in R]

    Rb, Zb, Lb, dψdn_R = fixed_boundary(shot, bnd)

    # ψ from image and coil currents
    Threads.@threads for i in eachindex(R)
        Vb = zero(Lb)
        @inbounds r = R[i]
        for j in eachindex(Z)
            @inbounds z = Z[j]
            # subtract image ψ
            Vb .= dψdn_R .* Green.(Rb, Zb, r, z)
            @inbounds ψ_f2f[j, i] -= -trapz(Lb, Vb)
            # add coil ψ
            @inbounds ψ_f2f[j, i] += μ₀ * Bp_fac * sum(coil.current * Green(coil, r, z) for coil in coils)
        end
    end

    return ψ_f2f
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

function chunk_indices(dims::Tuple{Vararg{Int}}, N::Int)
    # Total number of elements
    total_elements = prod(dims)

    # Calculate chunk size
    chunk_size, remainder = divrem(total_elements, N)

    # Split indices into N chunks
    chunks = []
    start_idx = 1
    for i in 1:N
        end_idx = start_idx + chunk_size - 1
        # Distribute the remainder among the first few chunks
        if i <= remainder
            end_idx += 1
        end
        chunk_1d = start_idx:end_idx
        chunk_multi = (CartesianIndices(dims)[i] for i in chunk_1d)
        push!(chunks, chunk_multi)
        start_idx = end_idx + 1
    end

    return chunks
end

# ******************************************
# Plots to check solution
# ******************************************
function check_fixed_eq_currents(
    EQfixed,
    coils::AbstractVector{<:AbstractCoil{T,C}},
    EQfree::Union{MXHEquilibrium.AbstractEquilibrium,Nothing}=nothing;
    resolution=257,
    Rmin=nothing,
    Rmax=nothing,
    Zmin=nothing,
    Zmax=nothing) where {T<:Real,C<:Real}

    Rmin0, Rmax0, Zmin0, Zmax0 = bounds(coils)
    if Rmin === nothing
        Rmin = Rmin0
    end
    if Rmax === nothing
        Rmax = Rmax0
    end
    if Zmin === nothing
        Zmin = Zmin0
    end
    if Zmax === nothing
        Zmax = Zmax0
    end

    R = range(Rmin, Rmax; length=resolution)
    Z = range(Zmin, Zmax; length=resolution)

    # ψ from fixed-boundary gEQDSK
    # make ψ at boundary zero, and very small value outside for plotting
    ψ0_fix, ψb_fix = MXHEquilibrium.psi_limits(EQfixed)
    ψb_fix = MXHEquilibrium.psi_boundary(EQfixed; precision=0.0)
    if ψb_fix === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        ψb_fix = MXHEquilibrium.psi_boundary(EQfixed)
    end
    σ₀ = sign(ψ0_fix - ψb_fix)
    ψ_fix = [EQfixed(r, z) for z in Z, r in R] .- ψb_fix
    ψ_fix = ifelse.(σ₀ * ψ_fix .> 0, ψ_fix, 1e-6 * ψ_fix)

    # ψ at the boundary is determined by the value of the currents
    # calculated in fixed_eq_currents
    ψ_f2f = fixed2free(EQfixed, coils, R, Z)

    # scale for plotting
    # this may get shifted boundary flux stays at the midpoint
    ψmax = 2.0 * abs(ψb_fix - ψ0_fix) * 0.99
    lvls = collect(-ψmax:0.1*ψmax:ψmax)
    clim = (-ψmax, ψmax)

    # Plot
    if isnothing(EQfree)
        # Overlay contours
        p = contour(R, Z, ψ_fix; levels=lvls, aspect_ratio=:equal,
            clim=clim, colorbar=false,
            linewidth=3, linecolor=:black,
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))

        # subtract ψ_f2f value at ψ_fix boundary to lineup contours
        ψtmp = ifelse.(σ₀ * ψ_fix .> 0, ψ_f2f, 0.0)
        offset = (σ₀ > 0 ? minimum(ψtmp) : maximum(ψtmp))
        contour!(R, Z, ψ_f2f .- offset; levels=lvls, colorbar=false,
            linewidth=1, linecolor=:red)
    else
        # Heat maps for free, fix, fix->free, and difference
        # ψ from free-boundary gEQDSK
        ψb_free = psi_boundary(EQfree; precision=0.0)
        if ψb_free === nothing
            # if the original boundary specified in EQfixed does not close, then find LCFS boundary
            ψb_free = psi_boundary(EQfree)
        end
        ψ_free = [EQfree(r, z) for z in Z, r in R]
        offset = ψb_free - ψb_fix

        lvls_off = lvls .+ offset
        clim_off = (minimum(lvls_off), maximum(lvls_off))

        pfree = heatmap(R, Z, ψ_free; clim=clim_off, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Free Boundary", ylabel="Z (m)",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_free; levels=lvls_off, linecolor=:black)

        pfix = heatmap(R, Z, ψ_fix; clim=clim, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Fixed Boundary",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_fix; levels=lvls, linecolor=:black)

        pf2f = heatmap(R, Z, ψ_f2f; clim=clim_off, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Calculated", xlabel="R (m)", ylabel="Z (m)",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_f2f; levels=lvls_off, linecolor=:black)

        pdiff = heatmap(R, Z, ψ_f2f - ψ_free; clim=clim, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Difference", xlabel="R (m)",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))

        contour!(R, Z, ψ_f2f - ψ_free; levels=lvls, linecolor=:black)
        p = plot(pfree, pfix, pf2f, pdiff; layout=(2, 2), size=(500, 550))
    end
    return p
end

"""
    coils_flux(Bp_fac::Real, coils::AbstractVector{<:AbstractCoil{T,C}}, R::AbstractVector{T}, Z::AbstractVector{T}) where {T<:Real,C<:Real}

Calculate flux from coils on a R, Z grid
"""
function coils_flux(Bp_fac::Real, coils::AbstractVector{<:AbstractCoil{T,C}}, R::AbstractVector{T}, Z::AbstractVector{T}) where {T<:Real,C<:Real}
    ψ = zeros(length(R), length(Z))
    @threads for i in eachindex(R)
        @inbounds r = R[i]
        for j in eachindex(Z)
            @inbounds z = Z[j]
            @inbounds ψ[i, j] += μ₀ * Bp_fac * sum(coil.current * Green(coil, r, z) for coil in coils)
        end
    end
    return ψ
end

function plot_coil_flux(
    Bp_fac::Real,
    coils::AbstractVector{<:AbstractCoil{T,C}},
    ψbound=0.0;
    resolution=129,
    clim=nothing,
    Rmin=nothing,
    Rmax=nothing,
    Zmin=nothing,
    Zmax=nothing
) where {T<:Real,C<:Real}

    Rmin0, Rmax0, Zmin0, Zmax0 = bounds(coils)
    if Rmin === nothing
        Rmin = Rmin0
    end
    if Rmax === nothing
        Rmax = Rmax0
    end
    if Zmin === nothing
        Zmin = Zmin0
    end
    if Zmax === nothing
        Zmax = Zmax0
    end

    R = range(Rmin, Rmax; length=resolution)
    Z = range(Zmin, Zmax; length=resolution)

    ψ = coils_flux(Bp_fac, coils, R, Z)

    if clim === nothing
        ψmax = maximum(abs, ψ)
        clim = (ψbound - ψmax, ψbound + ψmax)
    end

    # Plot heat maps for ψ from coil currents
    p = heatmap(R, Z, transpose(ψ);
        clim=clim,
        c=:diverging,
        aspect_ratio=:equal,
        linecolor=:black,
        title="Coil Flux",
        xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
    contour!(R, Z, transpose(ψ); levels=ψbound * [0.99, 1.00, 1.01], linecolor=:black)

    return p
end

# ########### #
# Convex Hull #
# ########### #
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
