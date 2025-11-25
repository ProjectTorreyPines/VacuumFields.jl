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


# in-place versions
@inline ψ!(result::AbstractVector{T}, source, R::T, Z::T; kwargs...) where {T<:Real}  = _pfunc!(Green, result, source, R, Z; kwargs...)
@inline dψ_dR!(result::AbstractVector{T}, source, R::T, Z::T; kwargs...) where {T<:Real} = _pfunc!(dG_dR, result, source, R, Z; kwargs...)
@inline dψ_dZ!(result::AbstractVector{T}, source, R::T, Z::T; kwargs...) where {T<:Real} = _pfunc!(dG_dZ, result, source, R, Z; kwargs...)
@inline d2ψ_dZ2!(result::AbstractVector{T}, source, R::T, Z::T; kwargs...) where {T<:Real} = _pfunc!(d2G_dZ2, result, source, R, Z; kwargs...)

# Faster in-place versions without kwargs
@inline ψ!(result::AbstractVector{T}, source, R::T, Z::T, Bp_fac::T) where {T<:Real} = _pfunc!(Green, result, source, R, Z, Bp_fac)
@inline dψ_dR!(result::AbstractVector{T}, source, R::T, Z::T, Bp_fac::T) where {T<:Real} = _pfunc!(dG_dR, result, source, R, Z, Bp_fac)
@inline dψ_dZ!(result::AbstractVector{T}, source, R::T, Z::T, Bp_fac::T) where {T<:Real} = _pfunc!(dG_dZ, result, source, R, Z, Bp_fac)
@inline d2ψ_dZ2!(result::AbstractVector{T}, source, R::T, Z::T, Bp_fac::T) where {T<:Real} = _pfunc!(d2G_dZ2, result, source, R, Z, Bp_fac)

"""Loop over coils and write ψ values into pre-allocated result array"""
function _pfunc!(Gfunc, result::AbstractVector{T}, coils::AbstractVector{<:Union{AbstractSingleCoil, IMAScoil}}, R::Real, Z::Real;
            COCOS::MXHEquilibrium.COCOS=MXHcocos11, Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp) where {T<:Real}
    for (i, coil) in enumerate(coils)
        coil_current_per_turn = current_per_turn(coil)::T
        if coil_current_per_turn == zero(T) 
            result[i] = zero(T)
        else
            result[i] = μ₀ * Bp_fac * Gfunc(coil, R, Z) * coil_current_per_turn
        end
    end
    return result
end

"""Loop over coils and write ψ values into pre-allocated result array"""
# function _pfunc!(Gfunc, result::AbstractVector{T}, coils::AbstractVector{<:Union{AbstractSingleCoil, IMAScoil}}, R::Real, Z::Real, Bp_fac::Real) where {T<:Real}
function _pfunc!(Gfunc, result::AbstractVector{T}, coils::AbstractVector{<:Union{AbstractSingleCoil, IMAScoil}}, R::T, Z::T, Bp_fac::T) where {T<:Real}
    
    for (i, coil) in enumerate(coils)
        coil_current_per_turn = current_per_turn(coil)::T
        if coil_current_per_turn == zero(T) 
            result[i] = zero(T)
        else
            result[i] = μ₀ * Bp_fac * Gfunc(coil, R, Z) * coil_current_per_turn
        end
    end
    return result
end


"""Loop over MultiCoils and write ψ values into pre-allocated result array"""
@inline function _pfunc!(Gfunc::Function, result::AbstractVector, mcoils::AbstractVector{<:MultiCoil}, R::Real, Z::Real, Bp_fac::Real; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11))
    for (i, mcoil) in enumerate(mcoils)
        coil_current_per_turn = current_per_turn(mcoil)
        if coil_current_per_turn == 0.0
            result[i] = 0.0
        else
            result[i] = sum(_pfunc(Gfunc, coil, R, Z; COCOS, Bp_fac, coil_current_per_turn) * mcoil.orientation[k] for (k, coil) in enumerate(mcoil.coils))
        end
    end
    return result
end


"""Loop over MultiCoils and write ψ values into pre-allocated result array"""
@inline function _pfunc!(Gfunc::Function, result::AbstractVector, mcoils::AbstractVector{<:MultiCoil}, R::Real, Z::Real;
                        COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    for (i, mcoil) in enumerate(mcoils)
        coil_current_per_turn = current_per_turn(mcoil)
        if coil_current_per_turn == 0.0
            result[i] = 0.0
        else
            result[i] = sum(_pfunc(Gfunc, coil, R, Z; COCOS, Bp_fac, coil_current_per_turn) * mcoil.orientation[k] for (k, coil) in enumerate(mcoil.coils))
        end
    end
    return result
end

"""Loop over images and write ψ values into pre-allocated result array"""
@inline function _pfunc!(Gfunc::Function, result::AbstractVector, images::AbstractVector{<:Image}, R::Real, Z::Real)
    for (i, image) in enumerate(images)
        Vb = (k, xx) -> image.dψdn_R[k] * Gfunc(image.Rb[k], image.Zb[k], R, Z)
        result[i] = -trapz(image.Lb, Vb)
    end
    return result
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
