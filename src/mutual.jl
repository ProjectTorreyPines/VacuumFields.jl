# mutual inductance between coils
# N.B. The flux through a coil is -σᵦₚ * (2π) ^ (1 - ϵᵦₚ) times ψ [or -σᵦₚ times poloidal flux ]
#      For coils: this cancels out the COCOS dependence in flux, leaving a -2π factor
#      For image: need the whole factor since the COCOS is buried into the image flux

"""
    mutual(C1::Union{AbstractCoil, IMASelement}, C2::PointCoil)

Compute the mutual inductance between an arbitrary coil or `IMAS.pf_active__coil___element` and a PointCoil
"""
function mutual(C1::Union{AbstractCoil, IMASelement}, C2::PointCoil; kwargs...)
    fac = -2π * μ₀ * turns(C2)
    return fac * Green(C1, C2.R, C2.Z; kwargs...)
end

"""
    mutual(C1::Union{AbstractCoil, IMASelement}, C2::DistributedCoil)

Compute the mutual inductance between an arbitrary coil or `IMAS.pf_active__coil___element` and a DistributedCoil
"""
function mutual(C1::Union{AbstractCoil, IMASelement}, C2::DistributedCoil; kwargs...)
    fac = -2π * μ₀ * turns(C2)
    return fac * sum(Green(C1, C2.R[k], C2.Z[k]; kwargs...) for k in eachindex(C2.R)) / length(C2.R)
end

"""
    mutual(C1::Union{AbstractCoil, IMASelement}, C2::Union{ParallelogramCoil, QuadCoil, IMASelement}; xorder::Int=3, yorder::Int=3)

Compute the mutual inductance between an arbitrary coil or `IMAS.pf_active__coil___element` and
    a ParallelogramCoil, QuadCoil, or `IMAS.pf_active__coil___element`

`xorder` and `yorder` give the order of Gauss-Legendre quadrature for integration over the coil area
"""
function mutual(C1::Union{AbstractCoil, IMASelement}, C2::Union{ParallelogramCoil, QuadCoil, IMASelement}; xorder::Int=3, yorder::Int=3)
    fac = -2π * μ₀ * turns(C2)
    f = (r, z) -> Green(C1, r, z; xorder = xorder+1, yorder = yorder+1)
    return fac * integrate(f, C2; xorder, yorder) / area(C2)
end

function mutual(C1::Union{AbstractCoil, IMASelement}, mcoil::MultiCoil; kwargs...)
    M = mutual(C1, mcoil.coils[1]; kwargs...) * mcoil.orientation[1]
    for k in eachindex(mcoil.coils)[2:end]
        M += mutual(C1, mcoil.coils[k]; kwargs...) * mcoil.orientation[k]
    end
    return M
    #return sum(mutual(C1, C2; kwargs...) for C2 in mcoil.coils)
end

"""
    mutual(C1::Union{AbstractCoil, IMASelement}, C2::IMAScoil; xorder::Int=3, yorder::Int=3)

Compute the mutual inductance between an arbitrary coil or `IMAS.pf_active__coil___element` and a `IMAS.pf_active__coil`

`xorder` and `yorder` give the order of Gauss-Legendre quadrature for integration over the coil area
"""
function mutual(C1::Union{AbstractCoil, IMASelement}, C2::IMAScoil; xorder::Int=3, yorder::Int=3)
    return sum(mutual(C1, element; xorder, yorder) for element in elements(C2))
end

"""
    mutual(C1::IMAScoil, C2::Union{AbstractCoil, IMAScoil, IMASelement}; xorder::Int=3, yorder::Int=3)

Compute the mutual inductance between an `IMAS.pf_active__coil` and an arbitrary coil, `IMAS.pf_active__coil___element`, or a `IMAS.pf_active__coil`

`xorder` and `yorder` give the order of Gauss-Legendre quadrature for integration over the coil area
"""
function mutual(C1::IMAScoil, C2::Union{AbstractCoil, IMAScoil, IMASelement}; xorder::Int=3, yorder::Int=3)
    return sum(mutual(element, C2; xorder, yorder) for element in elements(C1))
end


# Circuit

function mutual(C1::Union{AbstractCoil, IMASelement}, SC2::SeriesCircuit; kwargs...)
    return sum(SC2.signs[k] * mutual(C1, C2; kwargs...) for (k, C2) in enumerate(SC2.coils))
end

function mutual(SC1::SeriesCircuit, C2::Union{AbstractCoil, IMASelement, SeriesCircuit}; kwargs...)
    return sum(SC1.signs[k] * mutual(C1, C2; kwargs...) for (k, C1) in enumerate(SC1.coils))
end



# Plasma-coil mutuals

# Shifting plasma up by δZ is the same as shifting the coil down by δZ
function _pfunc(Pfunc::F1, image::Image, C::PointCoil, δZ; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11)) where {F1 <: Function}
    fac = -COCOS.sigma_Bp * (2π)^(1 - COCOS.exp_Bp)
    return fac * turns(C) * Pfunc(image, C.R, C.Z - δZ)
end

function _pfunc(Pfunc::F1, image::Image, C::DistributedCoil, δZ; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11)) where {F1 <: Function}
    fac = -COCOS.sigma_Bp * (2π)^(1 - COCOS.exp_Bp)
    return fac * turns(C) * sum(Pfunc(image, C.R[k], C.Z[k] - δZ) for k in eachindex(C.R)) / length(C.R)
end

function _pfunc(Pfunc::F1, image::Image, coil::IMAScoil, δZ;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                xorder::Int=default_order, yorder::Int=default_order) where {F1 <: Function}
    return sum(_pfunc(Pfunc, image, element, δZ; COCOS, xorder, yorder) for element in elements(coil))
end

function _pfunc(Pfunc::F1, image::Image, mcoil::MultiCoil, δZ;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                xorder::Int=default_order, yorder::Int=default_order) where {F1 <: Function}
    return sum(_pfunc(Pfunc, image, coil, δZ; COCOS, xorder, yorder) * mcoil.orientation[k] for (k, coil) in enumerate(mcoil.coils))
end

function _pfunc(Pfunc::F1, image::Image, C::Union{ParallelogramCoil, QuadCoil, IMASelement}, δZ;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                xorder::Int=default_order, yorder::Int=default_order) where {F1 <: Function}
    fac = -COCOS.sigma_Bp * (2π)^(1 - COCOS.exp_Bp)
    f = (r, z) -> Pfunc(image, r, z - δZ)
    return fac * turns(C) * integrate(f, C; xorder, yorder) / area(C)
end

"""
    mutual(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
        COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)

Compute the mutual inductance between an equilibrium and a coil,
    where the equilibrium is shifted vertically by `δZ`
"""
function mutual(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)
    return mutual(Image(EQ), C, plasma_current(EQ), δZ; COCOS, kwargs...)
end

"""
    mutual(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
           COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)

Compute the mutual inductance between an equilibrium's image current and a coil,
    where the equilibrium is shifted vertically by `δZ`
"""
function mutual(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)
    Ψ = _pfunc(ψ, image, C, δZ; COCOS, kwargs...)

    # negative sign since image flux is opposite plasma flux
    return - Ψ / Ip
end

"""
    dM_dZ(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
          COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)

Compute the Z derivative of the mutual inductance between an equilibrium and a coil,
    where the equilibrium is shifted vertically by `δZ`
"""
function dM_dZ(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
               COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)
    return dM_dZ(Image(EQ), C, plasma_current(EQ), δZ; COCOS, kwargs...)
end

"""
    dM_dZ(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
          COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)

Compute the Z derivative of the mutual inductance between an equilibrium's image current and a coil,
    where the equilibrium is shifted vertically by `δZ`
"""
function dM_dZ(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
               COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)
    # dψ/d(δZ) = -dψ_dZ
    dΨ_dZ = -_pfunc(dψ_dZ, image, C, δZ; COCOS, kwargs...)

    # negative sign since image flux is opposite plasma flux
    return -dΨ_dZ / Ip
end

"""
    d2M_dZ2(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
          COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)

Compute the second Z derivative of the mutual inductance between an equilibrium and a coil,
    where the equilibrium is shifted vertically by `δZ`
"""
function d2M_dZ2(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
                 COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)
    return d2M_dZ2(Image(EQ), C, plasma_current(EQ), δZ; COCOS, kwargs...)
end

"""
    d2M_dZ2(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
          COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)

Compute the second Z derivative of the mutual inductance between an equilibrium's image current and a coil,
    where the equilibrium is shifted vertically by `δZ`
"""
function d2M_dZ2(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
                 COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)
    # d²ψ/d(δZ)² = d²ψ_dZ²
    d2Ψ_dZ2 = _pfunc(d2ψ_dZ2, image, C, δZ; COCOS, kwargs...)

    # negative sign since image flux is opposite plasma flux
    return -d2Ψ_dZ2 / Ip
end

"""
    stability_margin(EQ::MXHEquilibrium.AbstractEquilibrium, coils::Vector{<:Union{AbstractCoil, IMAScoil}}; kwargs...)

Compute the m_s inductive stability margin for a given equilibrium and coils.
    Should be greater than 0.15 for vertical stability

First introduced in A. Portone, Nucl. Fusion 45 (2005) 926–932. https://doi.org/10.1088/0029-5515/45/8/021
"""
function stability_margin(EQ::MXHEquilibrium.AbstractEquilibrium, coils::Vector{<:Union{AbstractCoil, IMAScoil}}; kwargs...)
    return stability_margin(Image(EQ), coils, MXHEquilibrium.plasma_current(EQ); kwargs...)
end

"""
    stability_margin(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)

Compute the m_s inductive stability margin for a given equilibrium's image & plasma current and coils.
    Should be greater than 0.15 for vertical stability

First introduced in A. Portone, Nucl. Fusion 45 (2005) 926–932. https://doi.org/10.1088/0029-5515/45/8/021
"""
function stability_margin(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)
    b = [Ip * dM_dZ(image, C, Ip) for C in coils]
    Ktmp = (coil_current_per_turn, coil) -> coil_current_per_turn == 0.0 ? coil_current_per_turn : coil_current_per_turn * d2M_dZ2(image, coil, Ip)
    K = Ip * sum(Ktmp(current_per_turn(C), C) for C in coils)
    M = zeros(length(coils), length(coils))
    @threads for j in eachindex(coils)
        for k in eachindex(coils)
            k < j && continue
            M[k, j] = mutual(coils[k], coils[j]; xorder=order, yorder=order)
            (j != k) && (M[j, k] = M[k, j])
        end
    end
    # (bᵀ M⁻¹ b) / K - 1
    return dot(b, M \ b) / K - 1.0
end

"""
    normalized_growth_rate(EQ::MXHEquilibrium.AbstractEquilibrium, coils::Vector{<:Union{AbstractCoil, IMAScoil}}; kwargs...)

Compute the vertical growth rate (γ) and effective vertical time constant (τ, weighted L/R time) for a given equilibrium and coils

Return (γ, τ, γ * τ), where γ * τ < 10 for stability or controllability

This is the massless approximation and only use the passive conductors for computing τ (per advice from Olofsson)
"""
function normalized_growth_rate(EQ::MXHEquilibrium.AbstractEquilibrium, coils::Vector{<:Union{AbstractCoil, IMAScoil}}; kwargs...)
    return normalized_growth_rate(Image(EQ), coils, MXHEquilibrium.plasma_current(EQ); kwargs...)
end

"""
    normalized_growth_rate(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)

Compute the vertical growth rate (γ) and effective vertical time constant (τ, weighted L/R time), for a given
equilibrium's image & plasma current and coils

Return (γ, τ, γ * τ), where γ * τ < 10 for stability or controllability

This is the massless approximation and only use the passive conductors for computing τ (per advice from Olofsson)
"""
function normalized_growth_rate(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)
    b = [Ip * dM_dZ(image, C, Ip) for C in coils]
    Ktmp = (coil_current_per_turn, coil) -> coil_current_per_turn == 0.0 ? coil_current_per_turn : coil_current_per_turn * d2M_dZ2(image, coil, Ip)
    K = Ip * sum(Ktmp(current_per_turn(C), C) for C in coils)
    M = zeros(length(coils), length(coils))
    @threads for j in eachindex(coils)
        for k in eachindex(coils)
            k < j && continue
            M[k, j] = mutual(coils[k], coils[j]; xorder=order, yorder=order)
            (j != k) && (M[j, k] = M[k, j])
        end
    end

    Mstar = deepcopy(M)
    @threads for j in eachindex(coils)
        for k in eachindex(coils)[j:end]
            Mstar[k, j] -= b[k] * b[j] / K
            (j != k) && (Mstar[j, k] = Mstar[k, j])
        end
    end

    # reuse b vector for resistances
    @threads for j in eachindex(coils)
        b[j] = resistance(coils[j])
    end
    R = Diagonal(b)

    A0 = .- Mstar \ R

    D, V = eigen(A0)
    dmax = argmax(real.(D))
    γ = real(D[dmax])
    passive = [current_per_turn(C) == 0.0 for C in coils]
    @views v = V[passive, dmax]
    @views τ = dot(v, M[passive, passive], v) / dot(v, R[passive, passive], v)
    return γ, τ, γ * τ
end


