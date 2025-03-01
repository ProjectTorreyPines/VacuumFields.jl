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
    return coil_current == 0.0 ? coil_current : μ₀ * Bp_fac * Gfunc(coil, R, Z) * coil_current
end

@inline function _pfunc(Gfunc, mcoil::MultiCoil, R::Real, Z::Real;
    COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    cG = (coil_current, coil) -> coil_current == 0.0 ? coil_current : coil_current * Gfunc(coil, R, Z) 
    return μ₀ * Bp_fac * sum(cG(current(coil), coil) for coil in mcoil.coils)
end

# image flux
@inline function _pfunc(Gfunc, image::Image, R::Real, Z::Real)
    tid = Threads.threadid()
    image._Vb[tid] .= image.dψdn_R .* Gfunc.(image.Rb, image.Zb, R, Z)
    return -trapz(image.Lb, image._Vb[tid])
end

# series flux
@inline ψ(series::SeriesCircuit, R::Real, Z::Real; kwargs...)       = sum(ψ(coil, R, Z; kwargs...) for coil in series.coils)
@inline dψ_dR(series::SeriesCircuit, R::Real, Z::Real; kwargs...)   = sum(dψ_dR(coil, R, Z; kwargs...) for coil in series.coils)
@inline dψ_dZ(series::SeriesCircuit, R::Real, Z::Real; kwargs...)   = sum(dψ_dZ(coil, R, Z; kwargs...) for coil in series.coils)
@inline d2ψ_dZ2(series::SeriesCircuit, R::Real, Z::Real; kwargs...) = sum(d2ψ_dZ2(coil, R, Z; kwargs...) for coil in series.coils)