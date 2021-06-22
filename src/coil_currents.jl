function Green(X::Real, Y::Real, R::Real, Z::Real)
    XR = X*R
    m = 4.0*XR/((X+R)^2 + (Y-Z)^2) # this is k^2
    Km, Em = ellipk(m), ellipe(m)
    return (2.0*Em - (2.0 - m)*Km)*sqrt(XR/m)/(2π)
end

function fixed_boundary(Gfixed,skip=1)
    Sb = flux_surface(Gfixed, 0.)
    Rb, Zb = Sb.r[1:skip:end], Sb.z[1:skip:end]
    return (Rb, Zb)
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

function fixed_eq_currents(Gfixed, coils, skip=1)

    # boundary of fixed gEQDSK (radius, height, cumulative length)
    Rb, Zb = fixed_boundary(Gfixed,skip)
    
    Lb = cumlength(Rb,Zb)
    Bpol = poloidal_Bfield.(Gfixed, Rb, Zb) # Bpol = (1/R)*dPsi/dn
    
    # Calculate ψ from image currents on boundary at surface p near boundary
    ψ0,_ = psi_limits(Gfixed)
    Sp = flux_surface(Gfixed,1e-3*ψ0)
    Rp, Zp = Sp.r[1:end-1], Sp.z[1:end-1]
    Np = length(Rp)
    ψp = [@inbounds -trapz(Lb,Bpol.*Green.(Rb, Zb, Rp[j], Zp[j])) for j in 1:Np]

    # Compute coil currents needed to recreate ψ from image currents
    # Should use full coil geometry instead of singular currents
    Rc, Zc = unzip(coils)
    Nc = length(Rc)
    Gcp = μ₀.*[@inbounds Green(Rc[j], Zc[j], Rp[i], Zp[i]) for i in 1:Np, j in 1:Nc]
    Ic0 = (Gcp \ ψp)

    minimize_currents = false
    if minimize_currents
        λ₀ = sum((Gcp*Ic0 .- ψp).^2)*sum(Gcp*Ic0.^2)
        println(λ₀)
        function cost(Ic, λ=1000.)
            return sum((Gcp*Ic .- ψp).^2) + λ₀*λ*maximum(Ic.^2)
        end
        #println(cost(Ic0), Ic0*1e-3)
        Ic = Optim.minimizer(optimize(cost,Ic0))
        #println(cost(Ic), Ic*1e-3)
    end
    return Ic0
end

#******************************************
# Plots to check solution
#******************************************
function check_fixed_eq_currents(Gfixed,coils,currents,Gfree,KE::Union{Base.RefValue,Nothing}=nothing; skip=1)

    # mesh
    R, Z = Gfixed.r, Gfixed.z
    Rmin, Rmax = minimum(R), maximum(R)
    Zmin, Zmax = minimum(Z), maximum(Z)

    Rc, Zc = unzip(coils)

    # ψ from fixed-boundary gEQDSK
    ψfix = [min(Gfixed.psi_rz(r,z), 0.0) for z in Z, r in R]

    ψ0, _ = psi_limits(Gfixed)
    ψmax = abs(ψ0)*0.999
    lvls = -2ψmax:0.2*ψmax:2ψmax
    
    # ψ from free-boundary gEQDSK 
    _, ψb = psi_limits(Gfree)
    ψfree = [Gfree.psi_rz(r,z) for z in Z, r in R] .- ψb

    # ψ from image currents
    # boundary of fixed gEQDSK is at ψ=0
    Rb, Zb = fixed_boundary(Gfixed,skip)
    Lb = cumlength(Rb,Zb)
    Bpol = poloidal_Bfield.(Gfixed, Rb, Zb) # Bpol = (1/R)*dPsi/dn

    ψm = -[@inbounds trapz(Lb, Bpol.*Green.(Rb, Zb, r, z)) for z in Z, r in R]
        
    # ψ from coil currents
    # faster to do loop; I don't know why
    ψc = μ₀*[@inbounds sum(currents.*Green.(Rc, Zc, r, z)) for z in Z, r in R]

    # Plot
    p1 = heatmap(R, Z, ψfree, clim=(-2ψmax,2ψmax), c=:diverging,
                aspect_ratio=:equal,linecolor=:black,
                title="Free Boundary", ylabel="Z (m)",
                xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    contour!(R, Z, ψfree, levels=lvls, linecolor=:black)

    p2 = heatmap(R, Z, ψfix, clim=(-2ψmax,2ψmax), c=:diverging,
                aspect_ratio=:equal,linecolor=:black,
                title="Fixed Boundary",
                xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    # ψfix=0 doesn't plot properly since it's all zero outside boundary
    # this takes care of that
    ψtmp = [Gfixed.psi_rz(r,z) for z in Z, r in R]
    ψtmp = ifelse.(ψtmp.>0,ψtmp*1e-6,ψtmp)
    contour!(R, Z, ψtmp, levels=lvls, linecolor=:black)

    p3 = heatmap(R, Z, ψfix - ψm + ψc, clim=(-2ψmax,2ψmax), c=:diverging,
                aspect_ratio=:equal,linecolor=:black,
                title="Calculated", xlabel="R (m)", ylabel="Z (m)",
                xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    contour!(R, Z, ψfix - ψm + ψc, levels=lvls, linecolor=:black)

    p4 = heatmap(R, Z, ψfix - ψm + ψc - ψfree, clim=(-2ψmax,2ψmax), c=:diverging,
                aspect_ratio=:equal,linecolor=:black,
                title="Difference", xlabel="R (m)",
                xlim=(Rmin,Rmax),ylim=(Zmin,Zmax))
    contour!(R, Z, ψfix - ψm + ψc - ψfree, levels=lvls, linecolor=:black)

    p  = plot(p1, p2, p3, p4, layout=(2,2), size=(500,550))
    display(p)
    return
end