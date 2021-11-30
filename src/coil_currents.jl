abstract type AbstractCoil end

mutable struct PointCoil <: AbstractCoil
    R::Real
    Z::Real
end

"""
    ParallelogramCoil <: AbstractCoil

Parallelogram coil with the R, Z, ΔR, ΔZ, θ₁, θ₂ formalism (as used by EFIT, for example)
Here θ₁ and θ₂ are the shear angles along the x- and y-axes, respectively, in degrees.
"""
mutable struct ParallelogramCoil <: AbstractCoil
    R::Real
    Z::Real
    ΔR::Real
    ΔZ::Real
    θ₁::Real
    θ₂::Real
    spacing::Real
end

function ParallelogramCoil(R::Real, Z::Real, ΔR::Real, ΔZ::Real, θ₁::Real, θ₂::Real; spacing::Real=0.01)
    ParallelogramCoil(R, Z, ΔR, ΔZ, θ₁, θ₂, spacing)
end

mutable struct DistributedCoil <: AbstractCoil
    R::AbstractVector
    Z::AbstractVector
end

@Memoize.memoize function DistributedCoil(Rc::Real, Zc::Real, ΔR::Real, ΔZ::Real, θ₁::Real, θ₂::Real; spacing::Real=0.01)
    if spacing === nothing
        dR = LinRange(-0.5 * ΔR, 0.5 * ΔR, 2)
        dZ = LinRange(-0.5 * ΔZ, 0.5 * ΔZ, 2)
    else
        dR = LinRange(-0.5 * ΔR, 0.5 * ΔR, Int(floor(1.0 + ΔR / spacing)))
        dZ = LinRange(-0.5 * ΔZ, 0.5 * ΔZ, Int(floor(1.0 + ΔZ / spacing)))
    end
    α₁ = tan(π * θ₁ / 180.0)
    α₂ = tan(π * (θ₂ + 90.0) / 180.0)

    Nr = length(dR)
    Nz = length(dZ)
    N = Nr * Nz
    R, Z = zeros(N), zeros(N)
    for i in 1:Nr
        for j in 1:Nz
            k = j + (i - 1) * Nz
            R[k] = Rc + dR[i] - α₂ * dZ[j]
            Z[k] = Zc + dZ[j] + α₁ * dR[i]
        end
    end
    return DistributedCoil(R, Z)
end

function DistributedCoil(C::ParallelogramCoil)
    return DistributedCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ₁, C.θ₂; spacing=C.spacing)
end

# ============== #
#   Convex hull  #
# ============== #
function convex_hull(C::ParallelogramCoil)
    C = DistributedCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ₁, C.θ₂; spacing=nothing)
    return [[r,z] for (r, z) in zip(C.R, C.Z)]
end

function convex_hull(C::DistributedCoil)
    pts = [[r,z] for (r, z) in zip(C.R, C.Z)]
    return convex_hull(pts)
end

# ======== #
#   Plot   #
# ======== #
function plot_coil(C::ParallelogramCoil)
    return plot_coil(DistributedCoil(C))
end

function plot_coil(C::DistributedCoil)
    hull = convex_hull(C)
    plot!(VPolygon(hull), fillcolor=:black, alpha=0.2)
end

function plot_coil(C::PointCoil)
    plot!(Singleton([C.R,C.Z]), markercolor=:black)
end

function plot_coils(Cs::AbstractVector{T}) where {T <: AbstractCoil}
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
function Green(C::ParallelogramCoil, R::Real, Z::Real)
    return Green(DistributedCoil(C), R, Z)
end

function Green(C::DistributedCoil, R::Real, Z::Real)
    return sum(Green(x, y, R, Z) for (x, y) in zip(C.R, C.Z)) / length(C.R)
end

function Green(C::PointCoil, R::Real, Z::Real)
    return Green(C.R, C.Z, R, Z)
end

function Green(X::Real, Y::Real, R::Real, Z::Real)
    XR = X * R
    m = 4.0 * XR / ((X + R)^2 + (Y - Z)^2) # this is k^2
    if true # Use our own `Real` version of the elliptic functions to allow for ForwardDiff to work (copied from SpecialFunctions)
        Km = ellipk(m)
        Em = ellipe(m)
    else
        Km = SpecialFunctions.ellipk(m)
        Em = SpecialFunctions.ellipe(m)
    end
    return inv2π * (2.0 * Em - (2.0 - m) * Km) * sqrt(XR / m)
end

# ======== #
#   Utils  #
# ======== #
function cumlength(R, Z)
    # Length along boundary
    N = length(R)
    L = zeros(N)
    for i in 2:N
        @inbounds L[i] = L[i - 1] + sqrt((R[i] - R[i - 1])^2 + (Z[i] - Z[i - 1])^2)
    end
    return L
end

function bounds(Cs::AbstractVector{T}) where {T <: AbstractCoil}
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
function fixed_boundary(EQfixed)
    Sb = boundary(EQfixed)
    Rb, Zb = Sb.r, Sb.z
    Lb = cumlength(Rb, Zb)

    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    Bpol = poloidal_Bfield.(EQfixed, Rb, Zb) 
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    σ_ρθφ  = EQfixed.cocos.sigma_rhotp
    dψdn_R = σ_ρθφ * Bp_fac * Bpol
    return (Rb, Zb, Lb, dψdn_R)
end


"""
    ψp_on_fixed_eq_boundary(EQfixed,
                            fixed_coils=[],
                            fixed_currents=[],
                            ψbound=0.0;
                            Rx=[], Zx=[],
                            fraction_inside=0.999)

Calculate ψ from image currents on boundary at surface p near boundary
"""
function ψp_on_fixed_eq_boundary(EQfixed,
                                 fixed_coils=[],
                                 fixed_currents=[],
                                 ψbound=0.0;
                                 Rx=[], Zx=[],
                                 fraction_inside=0.999)
    ψ0, ψb = psi_limits(EQfixed)
    ψb = psi_boundary(EQfixed)

    # ψp is the flux from the image currents on the plasma boundary
    # which is equal and opposite to the flux from the plasma current
    # ψ₀ = ψp + ψim
    # ψ₁ = ψp + ψcoil
    # ψcoil = (ψ₁ - ψ₀) + ψim
    Sp = flux_surface(EQfixed, fraction_inside * (ψb - ψ0) + ψ0)
    Rp, Zp = Sp.r[1:end - 1], Sp.z[1:end - 1]
    append!(Rp, Rx)
    append!(Zp, Zx)
    ψp = zeros(length(Rp))

    # this is the image current contribution to the control points
    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)
    @threads for i = 1:length(Rp)
        ψp[i] = -trapz(Lb, dψdn_R .* Green.(Rb, Zb, Rp[i], Zp[i]))
    end

    # this is to account that the control points are inside of the LCFS
    # applied only to the plasma points
    ψp[1:end - length(Rx)] .-= (1.0 - fraction_inside) * (ψb - ψ0)

    # add in desired boundary flux
    ψbound != 0.0 && (ψp .+= ψbound)

    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp

    # Compute flux from fixed coils and subtract from ψp to match
    # This works whether ψp is a constant or a vector
    ψfixed = zeros(length(Rp))
    @threads for i = 1:length(Rp)
        for j in 1:length(fixed_coils)
            @inbounds ψfixed[i] += μ₀ * Bp_fac * Green(fixed_coils[j], Rp[i], Zp[i]) * fixed_currents[j]
        end
    end
    ψp = ψp .- ψfixed

    return Bp_fac, ψp, Rp, Zp
end

function currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils;
                              weights=nothing,
                              λ_regularize=1E-16,
                              return_cost=false)

    # Compute coil currents needed to recreate ψp at points (Rp,Zp)
    # Build matrix relating coil Green's functions to boundary points
    tp = typeof(sum([sum(c.R + c.Z) for c in coils]))
    Gcp = zeros(tp, length(Rp), length(coils))
    @threads for j = 1:length(coils)
        for i = 1:length(Rp)
            @inbounds Gcp[i,j] = μ₀ * Bp_fac * Green(coils[j], Rp[i], Zp[i])
        end
    end

    # handle weights
    if weights !== nothing
        if isa(weights, AbstractVector) weights = Diagonal(weights) end
        Gcp = weights * Gcp
        ψp  = weights * ψp
    end

    # Least-squares solve for coil currents
    if λ_regularize > 0
        # Least-squares with regularization
        # https://www.youtube.com/watch?v=9BckeGN0sF0
        reg_solve(A, b, λ) = inv(A' * A + λ * I) * A' * b
        Ic0 = reg_solve(Gcp, ψp, λ_regularize)
    else
        Ic0 = Gcp \ ψp
    end

    if return_cost
        cost(Ic) = norm(Gcp * Ic .- ψp) / μ₀ / length(Rp)
        return Ic0, cost(Ic0)
    else
        return Ic0
    end
end

function fixed_eq_currents(EQfixed,
                           coils,
                           fixed_coils=[],
                           fixed_currents=[],
                           ψbound=0.0;
                           λ_regularize=1E-16,
                           return_cost=false)

    Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(EQfixed, fixed_coils, fixed_currents, ψbound)

    return currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils;
                                λ_regularize=λ_regularize,
                                return_cost=return_cost)
end


# ***************************************************
# Transform fixed-boundary ψ to free-boundary ψ
# ***************************************************
function fixed2free(EQfixed, coils, currents, R, Z)
    
    ψ0, ψb = psi_limits(EQfixed)
    ψb = psi_boundary(EQfixed)
    σ₀ = sign(ψ0 - ψb)
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp

    tp = typeof(sum([sum(sum(coils[k].R) + sum(coils[k].Z) + sum(currents[k])) for k in 1:length(coils)]))
    ψ_f2f = [tp(EQfixed(r, z)) for z in Z, r in R] .- ψb
    ψ_f2f = ifelse.(σ₀ * ψ_f2f .> 0, ψ_f2f, 0) .+ ψb

    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)

    # ψ from image and coil currents
    @threads for i in 1:length(R)
        r = R[i]
        for j in 1:length(Z)
            z = Z[j]
            # subtract image ψ
            @inbounds ψ_f2f[j,i] -= -trapz(Lb, dψdn_R .* Green.(Rb, Zb, r, z))
            # add coil ψ
            @inbounds ψ_f2f[j,i] += μ₀ * Bp_fac * sum(currents .* Green.(coils, r, z))
        end
    end

    return ψ_f2f
end

# ******************************************
# Plots to check solution
# ******************************************
function check_fixed_eq_currents(EQfixed, coils, currents,
                                 EQfree::Union{AbstractEquilibrium,Nothing}=nothing;
                                 resolution=257,
                                 Rmin=nothing, Rmax=nothing, Zmin=nothing, Zmax=nothing)

    Rmin0, Rmax0, Zmin0, Zmax0 = bounds(coils)
    if Rmin === nothing Rmin = Rmin0 end
    if Rmax === nothing Rmax = Rmax0 end
    if Zmin === nothing Zmin = Zmin0 end
    if Zmax === nothing Zmax = Zmax0 end

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
    ψ_f2f = fixed2free(EQfixed, coils, currents, R, Z)

    # scale for plotting
    # this may get shifted boundary flux stays at the midpoint
    ψmax = 2.0 * abs(ψb_fix - ψ0_fix) * 0.99
    lvls = collect(-ψmax:0.1 * ψmax:ψmax)
    clim = (-ψmax, ψmax)

    # Plot
    if isnothing(EQfree)
        # Overlay contours
        p = contour(R, Z, ψ_fix, levels=lvls, aspect_ratio=:equal,
                    clim=clim, colorbar=false,
                    linewidth=3, linecolor=:black,
                    xlim=(Rmin, Rmax),ylim=(Zmin, Zmax))

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
                        aspect_ratio=:equal,linecolor=:black,
                        title="Free Boundary", ylabel="Z (m)",
                        xlim=(Rmin, Rmax),ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_free, levels=lvls_off, linecolor=:black)


        pfix = heatmap(R, Z, ψ_fix, clim=clim, c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Fixed Boundary",
                    xlim=(Rmin, Rmax),ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_fix, levels=lvls, linecolor=:black)

        pf2f = heatmap(R, Z, ψ_f2f, clim=clim_off, c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Calculated", xlabel="R (m)", ylabel="Z (m)",
                    xlim=(Rmin, Rmax),ylim=(Zmin, Zmax))
        contour!(R, Z, ψ_f2f, levels=lvls_off, linecolor=:black)

        pdiff = heatmap(R, Z, ψ_f2f - ψ_free, clim=clim, c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Difference", xlabel="R (m)",
                    xlim=(Rmin, Rmax),ylim=(Zmin, Zmax))

        contour!(R, Z, ψ_f2f - ψ_free, levels=lvls, linecolor=:black)
        p  = plot(pfree, pfix, pf2f, pdiff, layout=(2, 2), size=(500, 550))
    end
    return p
end

function plot_coil_flux(Bp_fac, coils, currents, ψbound=0.0;
                        resolution=257, clim=nothing,
                        Rmin=nothing, Rmax=nothing, Zmin=nothing, Zmax=nothing)

    Rmin0, Rmax0, Zmin0, Zmax0 = bounds(coils)
    if Rmin === nothing Rmin = Rmin0 end
    if Rmax === nothing Rmax = Rmax0 end
    if Zmin === nothing Zmin = Zmin0 end
    if Zmax === nothing Zmax = Zmax0 end

    R = range(Rmin, Rmax, length=resolution)
    Z = range(Zmin, Zmax, length=resolution)

    # ψ coil currents
    ψ = zeros(resolution, resolution)
    @threads for i in 1:resolution
        r = R[i]
        for j in 1:resolution
            z = Z[j]
            @inbounds ψ[j,i] += μ₀ * Bp_fac * sum(currents .* Green.(coils, r, z))
        end
    end

    if clim === nothing
        ψmax = maximum(abs.(ψ))
        clim = (ψbound - ψmax, ψbound + ψmax)
    end
    # Plot

    # Heat maps for ψ from coil currents
    p = heatmap(R, Z, ψ, clim=clim, c=:diverging,
                aspect_ratio=:equal,linecolor=:black,
                title="Coil Flux",
                xlim=(Rmin, Rmax),ylim=(Zmin, Zmax))
    contour!(R, Z, ψ, levels=ψbound * [0.99,1.00,1.01], linecolor=:black)

    return p
end
