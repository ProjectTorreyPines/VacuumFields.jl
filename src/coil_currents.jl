function Green(X::Real, Y::Real, R::Real, Z::Real)
    XR = X*R
    m = 4.0*XR/((X+R)^2 + (Y-Z)^2) # this is k^2
    Km, Em = ellipk(m), ellipe(m)
    return inv2π*(2.0*Em - (2.0 - m)*Km)*sqrt(XR/m)
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

function fixed_eq_currents(EQfixed, coils; minimize_currents=false)

    # Calculate ψ from image currents on boundary at surface p near boundary
    ψ0, ψb = psi_limits(EQfixed)
    Sp = flux_surface(EQfixed, 0.999*(ψb-ψ0) + ψ0)
    Rp, Zp = Sp.r[1:end-1], Sp.z[1:end-1]
    Np = length(Rp)
    ψp = zeros(Np)

    Rb, Zb, Lb, dψdn_R = fixed_boundary(EQfixed)
    @threads for i=1:Np
        ψp[i] = -trapz(Lb, dψdn_R .* Green.(Rb, Zb, Rp[i], Zp[i]))
    end

    # Compute coil currents needed to recreate ψ from image currents
    # Should use full coil geometry instead of singular currents
    Rc, Zc = unzip(coils)
    Nc = length(Rc)
    Gcp = zeros(Np, Nc)
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    @threads for j = 1:Nc
        for i=1:Np
            Gcp[i,j] = μ₀ * Bp_fac * Green(Rc[j], Zc[j], Rp[i], Zp[i])
        end
    end
    Ic0 = (Gcp \ ψp)

    if minimize_currents
        #λ₀ = sum((Gcp*Ic0 .- ψp).^2)*sum(Gcp*Ic0.^2)
        function cost(Ic, λ=0.001)
            return sum((Gcp*Ic .- ψp).^2) + λ*sum((μ₀*Ic).^2)
        end
        Ic = Optim.minimizer(optimize(cost,Ic0))
    else
        Ic = Ic0
    end

    return Ic
end

#***************************************************
# Transform fixed-boundary ψ to free-boundary ψ
#
#***************************************************
function fixed2free(EQfixed, coils, currents, R, Z)
    
    Bp_fac = EQfixed.cocos.sigma_Bp * (2π)^EQfixed.cocos.exp_Bp
    Rc, Zc = unzip(coils)

    ψ0, ψb = psi_limits(EQfixed)
    σ₀ = sign(ψ0-ψb)
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
            @inbounds ψ_f2f[j,i] += μ₀*Bp_fac*sum(currents.*Green.(Rc, Zc, r, z))
        end
    end

    return ψ_f2f
end

#******************************************
# Plots to check solution
#******************************************
function check_fixed_eq_currents(EQfixed, coils, currents,
                                 EQfree::Union{AbstractEquilibrium,Nothing}=nothing;
                                 resolution=257)

    Rc, Zc = unzip(coils)

    Rmin, Rmax = minimum(Rc), maximum(Rc)
    Zmin, Zmax = minimum(Zc), maximum(Zc)
    R = range(Rmin,Rmax,length=resolution)
    Z = range(Zmin,Zmax,length=resolution)

    # ψ from fixed-boundary gEQDSK
    ψ0_fix, ψb_fix = psi_limits(EQfixed)
    σ₀ = sign(ψ0_fix-ψb_fix)
    ψ_fix = [EQfixed(r,z) for z in Z, r in R] .- ψb_fix
    ψ_fix = ifelse.(σ₀*ψ_fix.>0, ψ_fix, 0) .+ ψb_fix
    
    ψ_f2f = fixed2free(EQfixed,coils,currents,R,Z)

    ψmax = 2.0*abs(ψb_fix - ψ0_fix)*0.99
    lvls = collect(-ψmax:0.1*ψmax:ψmax) .+ ψb_fix

    pfix = heatmap(R, Z, ψ_fix, clim=(-ψmax,ψmax), c=:diverging,
                   aspect_ratio=:equal,linecolor=:black,
                   title="Fixed Boundary",
                   xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    # ψ_fix=0 doesn't plot properly since it's all zero outside boundary
    # this takes care of that
    σ_Bp = EQfixed.cocos.sigma_Bp
    ψtmp = [EQfixed(r,z) for z in Z, r in R] .- ψb_fix
    ψtmp = ifelse.(σ₀*ψtmp.> 0, ψtmp, 1e-6*ψtmp) .+ ψb_fix
    contour!(R, Z, ψtmp, levels=lvls, linecolor=:black)

    pf2f = heatmap(R, Z, ψ_f2f, clim=(-ψmax,ψmax), c=:diverging,
                   aspect_ratio=:equal,linecolor=:black,
                   title="Calculated", xlabel="R (m)", ylabel="Z (m)",
                   xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    contour!(R, Z, ψ_f2f, levels=lvls, linecolor=:black)

    if isnothing(EQfree)
        p  = plot(pfix, pf2f, layout=(1,2), size=(550,250))
    else
        # ψ from free-boundary gEQDSK
         _, ψb_free = psi_limits(EQfree)
         ψ_free = [EQfree(r,z) for z in Z, r in R] .- ψb_free .+ ψb_fix

        # Plot
        pfree = heatmap(R, Z, ψ_free, clim=(-ψmax,ψmax), c=:diverging,
                        aspect_ratio=:equal,linecolor=:black,
                        title="Free Boundary", ylabel="Z (m)",
                        xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
        contour!(R, Z, ψ_free, levels=lvls, linecolor=:black)

        pdiff = heatmap(R, Z, ψ_f2f - ψ_free, clim=(-ψmax,ψmax), c=:diverging,
                    aspect_ratio=:equal,linecolor=:black,
                    title="Difference", xlabel="R (m)",
                    xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
        contour!(R, Z, ψ_f2f - ψ_free, levels=lvls, linecolor=:black)
        p  = plot(pfree, pfix, pf2f, pdiff, layout=(2,2), size=(500,550))
    end

    display(p)
    return 0.
end
