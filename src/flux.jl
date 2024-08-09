# trapezoidal rule
@inline trapz(x, y) = 0.5 * sum((x[k] - x[k-1]) * (y[k] + y[k-1]) for k in eachindex(x)[2:end])

# generalized functions
@inline ψ(source, R::Real, Z::Real; kwargs...)       = _pfunc(Green, source, R, Z; kwargs...)
@inline dψ_dR(source, R::Real, Z::Real; kwargs...)   = _pfunc(dG_dR, source, R, Z; kwargs...)
@inline dψ_dZ(source, R::Real, Z::Real; kwargs...)   = _pfunc(dG_dZ, source, R, Z; kwargs...)
@inline d2ψ_dZ2(source, R::Real, Z::Real; kwargs...) = _pfunc(d2G_dZ2, source, R, Z; kwargs...)


# coil flux
@inline function _pfunc(Gfunc, coil::Union{AbstractSingleCoil, IMAScoil}, R::Real, Z::Real;
                        COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp,
                        coil_current::Real=current(coil))
    return μ₀ * Bp_fac * Gfunc(coil, R, Z) * coil_current
end

@inline function _pfunc(Gfunc, mcoil::MultiCoil, R::Real, Z::Real;
    COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    return μ₀ * Bp_fac * sum(Gfunc(coil, R, Z) * current(coil) for coil in mcoil.coils)
end

# image flux
@inline function _pfunc(Gfunc, image::Image, R::Real, Z::Real)
    tid = Threads.threadid()
    image._Vb[tid] .= image.dψdn_R .* Gfunc.(image.Rb, image.Zb, R, Z)
    return -trapz(image.Lb, image._Vb[tid])
end