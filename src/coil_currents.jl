struct Coil{T <: AbstractVector, Q <: Integer}
    R :: T
    Z :: T
    Nr :: Q
    Nz :: Q
end

function coil(Rc::T,Zc::T,ΔR::T,ΔZ::T,θ₁::T,θ₂::T; spacing=0.01) where {T<:Real}
    dR = LinRange(-0.5*ΔR, 0.5*ΔR, Int(floor(1.0 + ΔR/spacing)))
    dZ = LinRange(-0.5*ΔZ, 0.5*ΔZ, Int(floor(1.0 + ΔZ/spacing)))
    α₁ = tan(π*θ₁/180.0)
    α₂ = tan(π*(θ₂+90.0)/180.0)

    Nr = length(dR)
    Nz = length(dZ)
    N = Nr*Nz
    R, Z = zeros(N), zeros(N)
    for i in 1:Nr
        for j in 1:Nz
            k = j + (i-1)*Nz
            R[k] = Rc + dR[i] - α₂*dZ[j]
            Z[k] = Zc + dZ[j] + α₁*dR[i]
        end
    end
    return Coil(R, Z, Nr, Nz)
end

function endpoints(C::Coil)
    idx = [1, C.Nz, C.Nr*C.Nz, (C.Nr-1)*C.Nz + 1, 1]
    return C.R[idx], C.Z[idx]
end

function bounds(Cs::AbstractVector{T}) where {T<:Coil}
    R,Z = endpoints(Cs[1])
    Rmin = minimum(R)
    Rmax = maximum(R)
    Zmin = minimum(Z)
    Zmax = maximum(Z)
    for C in Cs
        R, Z = endpoints(C)
        Rmin = min(Rmin,minimum(R))
        Rmax = max(Rmax,maximum(R))
        Zmin = min(Zmin,minimum(Z))
        Zmax = max(Zmax,maximum(Z))
    end
    return Rmin,Rmax,Zmin,Zmax
end

function plot_coils(Cs::AbstractVector{T}) where {T<:Coil}
    p = plot(aspect_ratio=:equal,legend=false)
    for C in Cs
        R,Z = endpoints(C)
        plot!(R,Z,linecolor=:black)
    end
    Rmin, Rmax, Zmin, Zmax = bounds(Cs)
    plot!(xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    display(p)
end

function Green(C::Coil, R::Real, Z::Real)
    return sum(Green(x, y, R, Z) for (x,y) in zip(C.R,C.Z))/(C.Nr*C.Nz)
end

function Green(C::Tuple{T,T}, R::Real, Z::Real) where {T<:Real}
    X, Y = C
    return Green(X, Y, R, Z)
end

function Green(X::Real, Y::Real, R::Real, Z::Real)
    XR = X*R
    m = 4.0*XR/((X+R)^2 + (Y-Z)^2 + 1E-16) # this is k^2
    if true
        Km, Em = SpecialFunctions.ellipk(m), SpecialFunctions.ellipe(m)
        return inv2π*(2.0*Em - (2.0 - m)*Km)*sqrt(XR/m)
    else
        K,E = Elliptic.ellipke(m)
        g = 2.0*E - (2.0 - m)*K
        g *= sqrt(XR/m)/(2π)
    end
end

function cumlength(R,Z)
    # Length along boundary
    N = length(R)
    L = zeros(N)
    for i in 2:N
        @inbounds L[i] = L[i-1] + sqrt((R[i] - R[i-1])^2 + (Z[i] - Z[i-1])^2)
    end
    return L
end

function bounds(Cs::Vector{Tuple{T,T}}) where {T<:Real}
    Rmin, Zmin = Cs[1]
    Rmax, Zmax = Cs[1]
    for i in 2:length(Cs)
        R,Z = Cs[i]
        Rmin = min(Rmin,R)
        Rmax = max(Rmax,R)
        Zmin = min(Zmin,Z)
        Zmax = max(Zmax,Z)
    end
    return Rmin, Rmax, Zmin, Zmax
end

function fixed_boundary(EQfixed)
    _,ψb = psi_limits(EQfixed)
    Sb = flux_surface(EQfixed, ψb)
    Rb, Zb = Sb.r, Sb.z
    Lb = cumlength(Rb,Zb)

    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    Bpol = poloidal_Bfield.(EQfixed, Rb, Zb) 
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    σ_ρθφ  = EQfixed.cocos.sigma_rhotp
    dψdn_R = σ_ρθφ*Bp_fac*Bpol
    return (Rb, Zb, Lb, dψdn_R)
end

function ψp_on_fixed_eq_boundary(EQfixed, ψbound=0.0)
    # Calculate ψ from image currents on boundary at surface p near boundary
    ψ0, ψb = psi_limits(EQfixed)
    Sp = flux_surface(EQfixed, 0.999*(ψb-ψ0) + ψ0)
    Rp, Zp = Sp.r[1:end-1], Sp.z[1:end-1]
    ψp = zeros(length(Rp))

    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)
    @threads for i=1:length(Rp)
        ψp[i] = -trapz(Lb, dψdn_R .* Green.(Rb, Zb, Rp[i], Zp[i]))
    end

    # add in desired boundary flux
    ψbound != 0.0 && (ψp .+= ψbound)

    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp

    return Bp_fac, ψp, Rp, Zp, ψbound
end

function currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, ψbound, coils;
    λ_minimize=0.0, λ_zerosum=0.0, λ_d3d_innersum=0.0,
    λ_regularize=1E-16,
    return_cost=false,
    verbose=false)

    # Compute coil currents needed to recreate ψ from image currents
    # Build matrix relating coil Green's functions to boundary points
    tp = typeof(sum([c[1]+c[2] for c in coils]))
    Gcp = zeros(tp, length(ψp), length(coils))
    
    @threads for j = 1:length(coils)
        for i=1:length(ψp)
            Gcp[i,j] = μ₀ * Bp_fac * Green(coils[j], Rp[i], Zp[i])
        end
    end

    # Least-squares solve for coil currents
    if λ_regularize>0
        # Least-squares with regularization
        # https://www.youtube.com/watch?v=9BckeGN0sF0
        reg_solve(A,b,λ) = inv(A'*A+λ*I)*A'*b
        Ic0 = reg_solve(Gcp, ψp, λ_regularize)
    else
        Ic0 = (Gcp \ ψp)
    end

    # use the least-squares regularized error as normalization for λ's part of optimization cost
    normalization = norm(Gcp*Ic0 .- ψp)/μ₀/length(ψp)

    function cost(Ic)
        (norm(Gcp*Ic .- ψp)/μ₀/length(ψp)
            + (λ_minimize*(μ₀*norm(Ic)/length(Ic))
               + λ_zerosum*(μ₀*abs(sum(Ic))/length(Ic))
              #+ λ_d3d_innersum*(μ₀*abs(sum(Ic[1:5])+sum(Ic[10:14]) + Ic[8] + Ic[17])/12.0)
              ) * normalization)
    end

    # Optional optimization:
    #    Total amplitude minimization
    #    Currents sum to zero
    if (λ_minimize>0.0 || λ_zerosum>0.0) || λ_d3d_innersum>0.0
        res = Optim.optimize(cost, Ic0, Optim.Newton(); autodiff=:forward)
        if verbose println(res) end
        Ic0 = Optim.minimizer(res)
    end
    
    if return_cost
        return Ic0, cost(Ic0)
    else
        return Ic0
    end
end

function fixed_eq_currents(EQfixed, coils, ψbound=0.0;
                           λ_minimize=0.0, λ_zerosum=0.0, λ_d3d_innersum=0.0,
                           λ_regularize=1E-16,
                           return_cost=false,
                           verbose=false)

    Bp_fac, ψp, Rp, Zp, ψbound = ψp_on_fixed_eq_boundary(EQfixed, ψbound)

    return currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, ψbound, coils;
                                λ_minimize=λ_minimize, λ_zerosum=λ_zerosum, λ_d3d_innersum=λ_d3d_innersum,
                                λ_regularize=λ_regularize,
                                return_cost=return_cost, verbose=verbose)
end

#***************************************************
# Transform fixed-boundary ψ to free-boundary ψ
#***************************************************
function fixed2free(EQfixed, coils, currents, R, Z)
    
    ψ0, ψb = psi_limits(EQfixed)
    σ₀ = sign(ψ0-ψb)
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp

    ψ_f2f = [EQfixed(r,z) for z in Z, r in R] .- ψb
    ψ_f2f = ifelse.(σ₀*ψ_f2f.>0, ψ_f2f, 0) .+ ψb

    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)

    # ψ from image and coil currents
    N = length(R)
    @threads for i in 1:N
        r = R[i]
        for j in 1:N
            z = Z[j]
            # subtract image ψ
            @inbounds ψ_f2f[j,i] -= -trapz(Lb,dψdn_R.*Green.(Rb, Zb, r, z))
            # add coil ψ
            @inbounds ψ_f2f[j,i] += μ₀*Bp_fac*sum(currents.*Green.(coils, r, z))
        end
    end

    return ψ_f2f
end

#******************************************
# Plots to check solution
#******************************************
function check_fixed_eq_currents(EQfixed, coils, currents,
                                 EQfree::Union{AbstractEquilibrium,Nothing}=nothing;
                                 resolution=257,
                                 Rmin=nothing, Rmax=nothing, Zmin=nothing, Zmax=nothing)

    Rmin0, Rmax0, Zmin0, Zmax0 = bounds(coils)
    if Rmin === nothing Rmin = Rmin0 end
    if Rmax === nothing Rmax = Rmax0 end
    if Zmin === nothing Zmin = Zmin0 end
    if Zmax === nothing Zmax = Zmax0 end

    R = range(Rmin,Rmax,length=resolution)
    Z = range(Zmin,Zmax,length=resolution)

    # ψ from fixed-boundary gEQDSK
    # make ψ at boundary zero, and very small value outside for plotting
    ψ0_fix, ψb_fix = psi_limits(EQfixed)
    σ₀ = sign(ψ0_fix-ψb_fix)
    ψ_fix = [EQfixed(r,z) for z in Z, r in R] .- ψb_fix
    ψ_fix = ifelse.(σ₀*ψ_fix.>0, ψ_fix, 1e-6*ψ_fix)
    
    # ψ at the boundary is determined by the value of the currents
    # calculated in fixed_eq_currents
    ψ_f2f = fixed2free(EQfixed,coils,currents,R,Z)

    # scale for plotting
    # this may get shifted boundary flux stays at the midpoint
    ψmax = 2.0*abs(ψb_fix - ψ0_fix)*0.99
    lvls = collect(-ψmax:0.1*ψmax:ψmax)
    clim = (-ψmax,ψmax)

    # Plot
    if isnothing(EQfree)
        # Overlay contours
        p = contour(R, Z, ψ_fix, levels=lvls, aspect_ratio=:equal,
                    clim=clim, colorbar=false,
                    linewidth=3, linecolor=:black,
                    xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))

        # subtract ψ_f2f value at ψ_fix boundary to lineup contours
        ψtmp = ifelse.(σ₀*ψ_fix.>0, ψ_f2f, 0.0)
        offset = (σ₀>0 ? minimum(ψtmp) : maximum(ψtmp))
        contour!(R, Z, ψ_f2f .- offset, levels=lvls, colorbar=false,
                 linewidth=1, linecolor=:red)
    else
        # Heat maps for free, fix, fix->free, and difference
        # ψ from free-boundary gEQDSK
        _, ψb_free = psi_limits(EQfree)
        ψ_free = [EQfree(r,z) for z in Z, r in R]
        offset = ψb_free - ψb_fix

        lvls_off = lvls .+ offset
        clim_off = (minimum(lvls_off), maximum(lvls_off))

        pfree = heatmap(R, Z, ψ_free, clim=clim_off, c=:diverging,
                        aspect_ratio=:equal,linecolor=:black,
                        title="Free Boundary", ylabel="Z (m)",
                        xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
        contour!(R, Z, ψ_free, levels=lvls_off, linecolor=:black)


        pfix = heatmap(R, Z, ψ_fix, clim=clim, c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Fixed Boundary",
                    xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
        contour!(R, Z, ψ_fix, levels=lvls, linecolor=:black)

        pf2f = heatmap(R, Z, ψ_f2f, clim=clim_off, c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Calculated", xlabel="R (m)", ylabel="Z (m)",
                    xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
        contour!(R, Z, ψ_f2f, levels=lvls_off, linecolor=:black)

        pdiff = heatmap(R, Z, ψ_f2f - ψ_free, clim=clim, c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Difference", xlabel="R (m)",
                    xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))

        contour!(R, Z, ψ_f2f - ψ_free, levels=lvls, linecolor=:black)
        p  = plot(pfree, pfix, pf2f, pdiff, layout=(2,2), size=(500,550))
    end

    return p
end
