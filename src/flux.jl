# generalized functions
@inline ψ(source, R::Real, Z::Real; kwargs...)       = _pfunc(Green, source, R, Z; kwargs...)
@inline dψ_dR(source, R::Real, Z::Real; kwargs...)   = _pfunc(dG_dR, source, R, Z; kwargs...)
@inline dψ_dZ(source, R::Real, Z::Real; kwargs...)   = _pfunc(dG_dZ, source, R, Z; kwargs...)
@inline d2ψ_dZ2(source, R::Real, Z::Real; kwargs...) = _pfunc(d2G_dZ2, source, R, Z; kwargs...)


# coil flux
@inline function _pfunc(Gfunc, coil::Union{AbstractSingleCoil, IMAScoil}, R::Real, Z::Real;
                        COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp,
                        coil_current_per_turn::Real=current_per_turn(coil))
    return coil_current_per_turn == 0.0 ? coil_current_per_turn : μ₀ * Bp_fac * Gfunc(coil, R, Z) * coil_current_per_turn
end

@inline function _pfunc(Gfunc, mcoil::MultiCoil, R::Real, Z::Real;
                        COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp,
                        coil_current_per_turn::Real=current_per_turn(mcoil))
    if coil_current_per_turn == 0.0
        return coil_current_per_turn
    end
    return sum(_pfunc(Gfunc, coil, R, Z; COCOS, Bp_fac, coil_current_per_turn) * mcoil.orientation[k] for (k, coil) in enumerate(mcoil.coils))
end

# image flux
@inline function _pfunc(Gfunc, image::Image, R::Real, Z::Real)
    Vb = (k, xx) -> image.dψdn_R[k] * Gfunc(image.Rb[k], image.Zb[k], R, Z)
    return -trapz(image.Lb, Vb)
end

# series flux
@inline ψ(series::SeriesCircuit, R::Real, Z::Real; kwargs...)       = sum(ψ(coil, R, Z; kwargs...) for coil in series.coils)
@inline dψ_dR(series::SeriesCircuit, R::Real, Z::Real; kwargs...)   = sum(dψ_dR(coil, R, Z; kwargs...) for coil in series.coils)
@inline dψ_dZ(series::SeriesCircuit, R::Real, Z::Real; kwargs...)   = sum(dψ_dZ(coil, R, Z; kwargs...) for coil in series.coils)
@inline d2ψ_dZ2(series::SeriesCircuit, R::Real, Z::Real; kwargs...) = sum(d2ψ_dZ2(coil, R, Z; kwargs...) for coil in series.coils)

"""
    flux_on_grid(Gtable::Array{<:Real, 3}, Rs::AbstractVector{T1}, Zs::AbstractVector{T2},
                 coils::Vector{<:Union{AbstractCoil, IMAScoil, IMASloop, IMASelement}};
                 COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                 Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
Calculate the magnetic flux on a (`Rs`, `Zs`) grid with Green's function table `Gtable` for the given coils and grid points.
"""
function flux_on_grid(Gtable::Array{T1, 3}, Rs::AbstractVector{T2}, Zs::AbstractVector{T3},
                      coils::Vector{<:Union{AbstractCoil, IMAScoil, IMASloop, IMASelement}};
                      COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                      Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp) where {T1 <: Real, T2 <: Real, T3 <: Real}
    Nr, Nz, _ = size(Gtable)
    Ψ = zeros(promote_type(T1, T2, T3), Nr, Nz)
    flux_on_grid!(Ψ, Gtable, Rs, Zs, coils; COCOS, Bp_fac)
    return Ψ
end

function flux_on_grid!(Ψ::Matrix{T}, Gtable::Array{<:Real, 3}, Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real},
                       coils::Vector{<:Union{AbstractCoil, IMAScoil, IMASloop, IMASelement}};
                       COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                       Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp) where {T <: Real}
    fill!(Ψ, zero(T))
    fac = Bp_fac * μ₀
    for k in eachindex(coils)
        Icpt = VacuumFields.current_per_turn(coils[k])
        Icpt == 0.0 && continue
        coeff = fac * Icpt
        @tturbo for j in eachindex(Zs), i in eachindex(Rs)
            Ψ[i, j] += coeff * Gtable[i, j, k]
        end
    end
end
