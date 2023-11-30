# trapezoidal rule
@inline trapz(x, y) = 0.5 * sum((x[k] - x[k-1]) * (y[k] + y[k-1]) for k in eachindex(x)[2:end])

# coil flux
@inline function ψ(coil::AbstractCoil, R::Real, Z::Real; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    return μ₀ * Bp_fac * Green(coil, R, Z) * coil.current
end

@inline function dψ_dR(coil::AbstractCoil, R::Real, Z::Real; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    return μ₀ * Bp_fac * dG_dR(coil, R, Z) * coil.current
end

@inline function dψ_dZ(coil::AbstractCoil, R::Real, Z::Real; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    return μ₀ * Bp_fac * dG_dZ(coil, R, Z) * coil.current
end

# image flux

@inline function ψ(image::Image, R::Real, Z::Real)
    tid = Threads.threadid()
    image._Vb[tid] .= image.dψdn_R .* Green.(image.Rb, image.Zb, R, Z)
    return -trapz(image.Lb, image._Vb[tid])
end

@inline function dψ_dR(image::Image, R::Real, Z::Real)
    tid = Threads.threadid()
    image._Vb[tid] .= image.dψdn_R .* dG_dR.(image.Rb, image.Zb, R, Z)
    return -trapz(image.Lb, image._Vb[tid])
end

@inline function dψ_dZ(image::Image, R::Real, Z::Real)
    tid = Threads.threadid()
    image._Vb[tid] .= image.dψdn_R .* dG_dZ.(image.Rb, image.Zb, R, Z)
    return -trapz(image.Lb, image._Vb[tid])
end