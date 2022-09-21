abstract type AbstractCoil end

mutable struct PointCoil{R1<:Real, R2<:Real, R3<:Real} <: AbstractCoil
    R::R1
    Z::R2
    current::R3
end

PointCoil(R, Z) = PointCoil(R, Z, 0.0)

"""
    ParallelogramCoil <: AbstractCoil

Parallelogram coil with the R, Z, ΔR, ΔZ, θ₁, θ₂ formalism (as used by EFIT, for example)
Here θ₁ and θ₂ are the shear angles along the x- and y-axes, respectively, in degrees.
"""
mutable struct ParallelogramCoil{R1<:Real, R2<:Real, R3<:Real, R4<:Real, R5<:Real, R6<:Real, UNR<:Union{Nothing,Real}, R7<:Real} <: AbstractCoil
    R::R1
    Z::R2
    ΔR::R3
    ΔZ::R4
    θ₁::R5
    θ₂::R6
    spacing::UNR
    current::R7
end

ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing) = ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing, 0.0)

function ParallelogramCoil(R::Real, Z::Real, ΔR::Real, ΔZ::Real, θ₁::Real, θ₂::Real; spacing::Union{Nothing,Real}=0.01)
    return ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing)
end

mutable struct DistributedCoil{AVR1<:AbstractVector{<:Real}, AVR2<:AbstractVector{<:Real}, R1<:Real} <: AbstractCoil
    R::AVR1
    Z::AVR2
    current::R1
end

DistributedCoil(R, Z) = DistributedCoil(R, Z, 0.0)

function DistributedParallelogramCoil(Rc::Real, Zc::Real, ΔR::Real, ΔZ::Real, θ₁::Real, θ₂::Real; spacing::Union{Nothing,Real}=0.01)
    if spacing === nothing
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
    return DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ₁, C.θ₂; spacing=C.spacing)
end

# ============== #
#   Convex hull  #
# ============== #
function convex_hull(C::ParallelogramCoil)
    C = DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ₁, C.θ₂; spacing=nothing)
    return [[r, z] for (r, z) in zip(C.R, C.Z)]
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
    plot!(Shape(R, Z), color=:black, alpha=0.2, label="")
end

function plot_coil(C::PointCoil)
    plot!([C.R], [C.Z], marker=:circle, markercolor=:black)
end

function plot_coils(Cs::AbstractVector{<:AbstractCoil})
    p = plot(aspect_ratio=:equal, legend=false)
    for C in Cs
        plot_coil(C)
    end
    Rmin, Rmax, Zmin, Zmax = bounds(Cs)
    Rmin -= 0.05 * (abs(Rmin) + abs(Rmax))
    Rmax += 0.05 * (abs(Rmin) + abs(Rmax))
    Zmin -= 0.05 * (abs(Zmin) + abs(Zmax))
    Zmax += 0.05 * (abs(Zmin) + abs(Zmax))
    plot!(xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
    display(p)
end

# ======== #
#   Green  #
# ======== #
function Green(C::ParallelogramCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return Green(DistributedCoil(C), R, Z, scale_factor)
end

function Green(C::DistributedCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return sum(Green(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end

function Green(C::PointCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return Green(C.R, C.Z, R, Z, scale_factor)
end

function Green(X::Real, Y::Real, R::Real, Z::Real, scale_factor::Real=1.0)
    XR = X * R
    m = 4.0 * XR / ((X + R)^2 + (Y - Z)^2) # this is k^2
    if true # Use our own `Real` version of the elliptic functions to allow for ForwardDiff to work (copied from SpecialFunctions)
        Km = ellipk(m)
        Em = ellipe(m)
    else
        Km = SpecialFunctions.ellipk(m)
        Em = SpecialFunctions.ellipe(m)
    end
    return inv2π * (2.0 * Em - (2.0 - m) * Km) * sqrt(XR / m) * scale_factor
end

# ======== #
#   Utils  #
# ======== #
function cumlength(R, Z)
    # Length along boundary
    L = Array{Float64,1}(undef, length(R))
    @inbounds L[1] = 0
    for i in 2:length(R)
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
function fixed_boundary(EQfixed::Equilibrium.AbstractEquilibrium)
    Sb = boundary(EQfixed)
    Rb, Zb = Sb.r, Sb.z
    Lb = cumlength(Rb, Zb)

    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    Bpol = poloidal_Bfield.(EQfixed, Rb, Zb)
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    σ_ρθφ = EQfixed.cocos.sigma_rhotp
    dψdn_R = σ_ρθφ * Bp_fac * Bpol
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
    EQfixed::Equilibrium.AbstractEquilibrium,
    fixed_coils::AbstractVector{<:AbstractCoil}=AbstractCoil[],
    ψbound::Real=0.0;
    Rx::AbstractVector{<:Real}=Real[],
    Zx::AbstractVector{<:Real}=Real[],
    fraction_inside::Union{Nothing,Real}=0.999)

    ψ0, ψb = psi_limits(EQfixed)
    ψb = psi_boundary(EQfixed)

    # ψp is the flux from the image currents on the plasma boundary
    # which is equal and opposite to the flux from the plasma current
    # ψ₀ = ψp + ψim
    # ψ₁ = ψp + ψcoil
    # ψcoil = (ψ₁ - ψ₀) + ψim
    Sp = flux_surface(EQfixed, fraction_inside * (ψb - ψ0) + ψ0)
    n = Int(floor(length(Sp.r) / 100.0)) + 1 # roughly at most 100 points
    Rp, Zp = Sp.r[1:n:end-1], Sp.z[1:n:end-1]

    append!(Rp, Rx)
    append!(Zp, Zx)

    # this is the image current contribution to the control points
    ψp = Array{Float64,1}(undef, length(Rp))
    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)
    @threads for i = 1:length(Rp)
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
    @threads for j in 1:length(fixed_coils)
        for i = 1:length(Rp)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rp[i], Zp[i]) * fixed_coils[j].current
        end
    end
    ψp = ψp .- ψfixed

    return Bp_fac, ψp, Rp, Zp
end

"""
    field_null_on_boundary(ψp_constant, Rp, Zp,
                           fixed_coils=[],
                           ψbound=0.0,
                           cocos=11)

Account for effect of fixed coils on ψp_constant
"""
function field_null_on_boundary(ψp_constant::Real,
    Rp::AbstractVector,
    Zp::AbstractVector,
    fixed_coils::AbstractVector{<:AbstractCoil}=AbstractCoil[],
    ψbound::Real=0.0,
    cocos::Int=11)

    # add in desired boundary flux
    ψbound != 0.0 && (ψp_constant .+= ψbound)

    Bp_fac = Equilibrium.cocos(cocos).sigma_Bp * (2π)^Equilibrium.cocos(cocos).exp_Bp

    # Compute flux from fixed coils and subtract from ψp to match
    ψfixed = zeros(length(Rp))
    @threads for j in 1:length(fixed_coils)
        for i = 1:length(Rp)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rp[i], Zp[i]) * fixed_coils[j].current
        end
    end
    ψp = ψp_constant .- ψfixed

    return Bp_fac, ψp, Rp, Zp
end

function currents_to_match_ψp(
    Bp_fac::Float64,
    ψp::Vector{<:Real},
    Rp::Vector{<:Real},
    Zp::Vector{<:Real},
    coils::AbstractVector{<:AbstractCoil};
    weights::Vector{Float64}=Float64[],
    λ_regularize::Float64=1E-16,
    return_cost::Bool=false,
    tp=Float64)

    # Compute coil currents needed to recreate ψp at points (Rp,Zp)
    # Build matrix relating coil Green's functions to boundary points
    Gcp = Array{tp,2}(undef, length(Rp), length(coils))
    @threads for j = 1:length(coils)
        for i = 1:length(Rp)
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
    for k in 1:length(coils)
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
    EQfixed::Equilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:AbstractCoil},
    fixed_coils::AbstractVector{<:AbstractCoil}=AbstractCoil[],
    ψbound::Real=0.0;
    λ_regularize::Real=1E-16,
    return_cost::Bool=false)

    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(EQfixed, fixed_coils, ψbound)
    return currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; λ_regularize, return_cost)
end


# ***************************************************
# Transform fixed-boundary ψ to free-boundary ψ
# ***************************************************
"""
    fixed2free(
        EQfixed::Equilibrium.AbstractEquilibrium,
        n_point_coils::Integer,
        R::AbstractVector,
        Z::AbstractVector)

Distribute n point coils around fixed boundary plasma to get a free boundary ψ map
"""
function fixed2free(
    EQfixed::Equilibrium.AbstractEquilibrium,
    n_coils::Integer;
    Rx::AbstractVector{<:Real}=Real[],
    Zx::AbstractVector{<:Real}=Real[],
    fraction_inside::Union{Nothing,Real}=0.999,
    Rgrid::AbstractVector{<:Real}=EQfixed.r,
    Zgrid::AbstractVector{<:Real}=EQfixed.z)

    coils = encircling_coils(EQfixed, n_coils)
    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(EQfixed, coils; fraction_inside, Rx, Zx)
    currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; λ_regularize=1E-15)
    return transpose(fixed2free(EQfixed, coils, Rgrid, Zgrid))
end

function encircling_coils(EQfixed::Equilibrium.AbstractEquilibrium, n_coils::Integer)
    bnd = Equilibrium.boundary(EQfixed)
    R0 = sum(bnd.r) / length(bnd.r)
    Z0 = sum(bnd.z) / length(bnd.z)
    t = LinRange(0, 2π, n_coils + 1)[1:n_coils]
    a = R0 * 0.99
    b = (maximum(bnd.z) - minimum(bnd.z)) * 1.5
    b = max(a, b)
    a = min(a, b)
    coils_r = a .* cos.(t) .+ R0
    coils_z = b .* sin.(t) .+ Z0
    return [PointCoil(r, z) for (r, z) in zip(coils_r, coils_z)]
end

function fixed2free(
    EQfixed::Equilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:AbstractCoil},
    R::AbstractVector{<:Real},
    Z::AbstractVector{<:Real};
    tp=Float64)

    ψb = psi_boundary(EQfixed)
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp

    Sb = boundary(EQfixed)
    ψ_f2f = [in_boundary(Sb, (r, z)) ? tp(EQfixed(r, z)) : ψb for z in Z, r in R]

    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)

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

function fixed2free(
    EQfixed::Equilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:ParallelogramCoil},
    R::AbstractVector{<:Real},
    Z::AbstractVector{<:Real};
    tp=Float64)
    dcoils = DistributedCoil.(coils)
    return fixed2free(EQfixed, dcoils, R, Z; tp)
end

# ******************************************
# Plots to check solution
# ******************************************
function check_fixed_eq_currents(
    EQfixed,
    coils::AbstractVector{<:AbstractCoil},
    EQfree::Union{AbstractEquilibrium,Nothing}=nothing;
    resolution=257,
    Rmin=nothing,
    Rmax=nothing,
    Zmin=nothing,
    Zmax=nothing)

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

    R = range(Rmin, Rmax, length=resolution)
    Z = range(Zmin, Zmax, length=resolution)

    # ψ from fixed-boundary gEQDSK
    # make ψ at boundary zero, and very small value outside for plotting
    ψ0_fix, ψb_fix = psi_limits(EQfixed)
    ψb_fix = psi_boundary(EQfixed)
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
        p = contour(R, Z, ψ_fix, levels=lvls, aspect_ratio=:equal,
            clim=clim, colorbar=false,
            linewidth=3, linecolor=:black,
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))

        # subtract ψ_f2f value at ψ_fix boundary to lineup contours
        ψtmp = ifelse.(σ₀ * ψ_fix .> 0, ψ_f2f, 0.0)
        offset = (σ₀ > 0 ? minimum(ψtmp) : maximum(ψtmp))
        contour!(R, Z, ψ_f2f .- offset, levels=lvls, colorbar=false,
            linewidth=1, linecolor=:red)
    else
        # Heat maps for free, fix, fix->free, and difference
        # ψ from free-boundary gEQDSK
        ψb_free = psi_boundary(EQfree)
        ψ_free = [EQfree(r, z) for z in Z, r in R]
        offset = ψb_free - ψb_fix

        lvls_off = lvls .+ offset
        clim_off = (minimum(lvls_off), maximum(lvls_off))

        pfree = heatmap(R, Z, ψ_free, clim=clim_off, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Free Boundary", ylabel="Z (m)",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_free, levels=lvls_off, linecolor=:black)

        pfix = heatmap(R, Z, ψ_fix, clim=clim, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Fixed Boundary",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_fix, levels=lvls, linecolor=:black)

        pf2f = heatmap(R, Z, ψ_f2f, clim=clim_off, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Calculated", xlabel="R (m)", ylabel="Z (m)",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_f2f, levels=lvls_off, linecolor=:black)

        pdiff = heatmap(R, Z, ψ_f2f - ψ_free, clim=clim, c=:diverging,
            aspect_ratio=:equal, linecolor=:black,
            title="Difference", xlabel="R (m)",
            xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))

        contour!(R, Z, ψ_f2f - ψ_free, levels=lvls, linecolor=:black)
        p = plot(pfree, pfix, pf2f, pdiff, layout=(2, 2), size=(500, 550))
    end
    return p
end

"""
    coils_flux(Bp_fac, coils, R, Z)

Calculate flux from coils on a R, Z grid
"""
function coils_flux(Bp_fac::Real, coils::AbstractVector{<:AbstractCoil}, R::AbstractVector{<:Real}, Z::AbstractVector{<:Real})
    ψ = zeros(length(R), length(Z))
    @threads for i in 1:length(R)
        @inbounds r = R[i]
        for j in 1:length(Z)
            @inbounds z = Z[j]
            @inbounds ψ[i, j] += μ₀ * Bp_fac * sum([coil.current for coil in coils] .* Green.(coils, r, z))
        end
    end
    return ψ
end

function plot_coil_flux(Bp_fac::Real, coils::AbstractVector{<:AbstractCoil}, ψbound=0.0;
    resolution=129, clim=nothing,
    Rmin=nothing, Rmax=nothing, Zmin=nothing, Zmax=nothing)

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

    R = range(Rmin, Rmax, length=resolution)
    Z = range(Zmin, Zmax, length=resolution)

    ψ = coils_flux(Bp_fac, coils, R, Z)

    if clim === nothing
        ψmax = maximum(abs.(ψ))
        clim = (ψbound - ψmax, ψbound + ψmax)
    end

    # Plot heat maps for ψ from coil currents
    p = heatmap(R, Z, transpose(ψ),
        clim=clim,
        c=:diverging,
        aspect_ratio=:equal,
        linecolor=:black,
        title="Coil Flux",
        xlim=(Rmin, Rmax), ylim=(Zmin, Zmax))
    contour!(R, Z, transpose(ψ), levels=ψbound * [0.99, 1.00, 1.01], linecolor=:black)

    return p
end

# ******************************************
# Convex Hull
# ******************************************
struct Point
    x::Float64
    y::Float64
end

function Base.isless(p::Point, q::Point)
    return p.x < q.x || (p.x == q.x && p.y < q.y)
end

function isrightturn(p::Point, q::Point, r::Point)::Bool
    return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x) < 0.0
end

function halfhull(points::Vector{Point})
    halfhull = points[1:2]
    for p in points[3:end]
        push!(halfhull, p)
        while length(halfhull) > 2 && !isrightturn(halfhull[end-2], halfhull[end-1], halfhull[end])
            deleteat!(halfhull, length(halfhull) - 1)
        end
    end
    return halfhull
end

function grahamscan(points::Vector{Point})
    sort!(points)
    upperhull = halfhull(points)
    lowerhull = halfhull(reverse!(points))
    return Point[upperhull; lowerhull[2:end-1]]
end

function convex_hull(xy_points::Vector{Point}; closed_polygon::Bool)
    tmp = [(k.x, k.y) for k in grahamscan(xy_points)]
    if closed_polygon
        return push!(tmp, tmp[1])
    else
        return tmp
    end
end

function convex_hull(xy::Vector; closed_polygon::Bool)
    xy_points = [Point(xx, yx) for (xx, yx) in xy]
    return convex_hull(xy_points; closed_polygon)
end

function convex_hull(x::Vector, y::Vector; closed_polygon::Bool)
    xy_points = [Point(xx, yx) for (xx, yx) in zip(x, y)]
    return convex_hull(xy_points; closed_polygon)
end
