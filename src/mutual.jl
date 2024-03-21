# mutual inductance between coils
# N.B. The flux through a coil is -σᵦₚ * (2π) ^ (1 - ϵᵦₚ) times ψ [or -σᵦₚ times poloidal flux ]
#      For coils: this cancels out the COCOS dependence in flux, leaving a -2π factor
#      For image: need the whole factor since the COCOS is buried into the image flux

function mutual(C1::Union{AbstractCoil, IMASelement}, C2::PointCoil)
    fac = -2π * μ₀ * turns(C1) * turns(C2)
    return fac * Green(C1, C2.R, C2.Z)
end

function mutual(C1::Union{AbstractCoil, IMASelement}, C2::DistributedCoil)
    fac = -2π * μ₀ * turns(C1) * turns(C2)
    return fac * sum(Green(C1, C2.R[k], C2.Z[k]) for k in eachindex(C2.R)) / length(C2.R)
end

function mutual(C1::Union{AbstractCoil, IMASelement}, C2::Union{ParallelogramCoil, QuadCoil, IMASelement}; xorder::Int=3, yorder::Int=3)
    fac = -2π * μ₀ * turns(C1) * turns(C2)
    f = (r, z) -> Green(C1, r, z; xorder = xorder+1, yorder = yorder+1)
    return fac * integrate(f, C2; xorder, yorder) / area(C2)
end

function mutual(C1::Union{AbstractCoil, IMASelement}, C2::IMAScoil; xorder::Int=3, yorder::Int=3)
    return sum(mutual(C1, element; xorder, yorder) for element in C2.element)
end

function mutual(C1::IMAScoil, C2::Union{AbstractCoil, IMAScoil, IMASelement}; xorder::Int=3, yorder::Int=3)
    return sum(mutual(element, C2; xorder, yorder) for element in C1.element)
end


# Plasma-coil mutuals

# Shifting plasma up by δZ is the same as shifting the coil down by δZ
function _pfunc(Pfunc, image::Image, C::PointCoil, δZ; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11))
    fac = -COCOS.sigma_Bp * (2π)^(1 - COCOS.exp_Bp)
    return fac * Pfunc(image, C.R, C.Z - δZ)
end

function _pfunc(Pfunc, image::Image, C::DistributedCoil, δZ; COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11))
    fac = -COCOS.sigma_Bp * (2π)^(1 - COCOS.exp_Bp)
    return fac * sum(Pfunc(image, C.R[k], C.Z[k] - δZ) for k in eachindex(C.R)) / length(C.R)
end

function _pfunc(Pfunc, image::Image, coil::IMAScoil, δZ;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                xorder::Int=default_order, yorder::Int=default_order)
    return sum(_pfunc(Pfunc, image, element, δZ; COCOS, xorder, yorder) for element in coil.element)
end

function _pfunc(Pfunc, image::Image, C::Union{ParallelogramCoil, QuadCoil, IMASelement}, δZ;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                xorder::Int=default_order, yorder::Int=default_order)
    fac = -COCOS.sigma_Bp * (2π)^(1 - COCOS.exp_Bp)
    f = (r, z) -> Pfunc(image, r, z - δZ)
    return fac * integrate(f, C; xorder, yorder) / area(C)
end

function mutual(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)
    return mutual(Image(EQ), C, plasma_current(EQ), δZ; COCOS, kwargs...)
end
function mutual(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)
    Ψ = _pfunc(ψ, image, C, δZ; COCOS, kwargs...)

    # negative sign since image flux is opposite plasma flux
    return -fac * turns(C) * Ψ / Ip
end

function dM_dZ(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
               COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)
    return dM_dZ(Image(EQ), C, plasma_current(EQ), δZ; COCOS, kwargs...)
end
function dM_dZ(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
               COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)
    # dψ/d(δZ) = -dψ_dZ
    dΨ_dZ = -_pfunc(dψ_dZ, image, C, δZ; COCOS, kwargs...)

    # negative sign since image flux is opposite plasma flux
    return -turns(C) * dΨ_dZ / Ip
end

function d2M_dZ2(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;
                 COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)
    return d2M_dZ2(Image(EQ), C, plasma_current(EQ), δZ; COCOS, kwargs...)
end
function d2M_dZ2(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;
                 COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)
    # d²ψ/d(δZ)² = d²ψ_dZ²
    d2Ψ_dZ2 = _pfunc(d2ψ_dZ2, image, C, δZ; COCOS, kwargs...)

    # negative sign since image flux is opposite plasma flux
    return -turns(C) * d2Ψ_dZ2 / Ip
end

# The m_s inductive stability margin
# (bᵀ M⁻¹ b) / K - 1
# Should be greater than 0.15
stability_margin(EQ, coils; kwargs...) = stability_margin(Image(EQ), coils, MXHEquilibrium.plasma_current(EQ); kwargs...)
function stability_margin(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)
    b = [Ip * dM_dZ(image, C, Ip) for C in coils]
    K = Ip * sum(current(C) * d2M_dZ2(image, C, Ip) for C in coils)
    M = zeros(length(coils), length(coils))
    for j in eachindex(coils)
        for k in eachindex(coils)
            k < j && continue
            M[k, j] = mutual(coils[k], coils[j]; xorder=order, yorder=order)
            (j != k) && (M[j, k] = M[k, j])
        end
    end
    return dot(b, M \ b) / K - 1.0
end


# Finds the vertical growth rate (γ) and effective vertical time constant (τ, weighted L/R time)
# γ * τ < 10 for stability
# Here we've implemented the massless approximation,
#    and only use the passive conductors for computing τ (per advice from Olofsson)
normalized_growth_rate(EQ, coils; kwargs...) = normalized_growth_rate(Image(EQ), coils, MXHEquilibrium.plasma_current(EQ); kwargs...)
function normalized_growth_rate(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)
    b = [Ip * dM_dZ(image, C, Ip) for C in coils]
    K = Ip * sum(current(C) * d2M_dZ2(image, C, Ip) for C in coils)
    M = zeros(length(coils), length(coils))
    for j in eachindex(coils)
        for k in eachindex(coils)
            k < j && continue
            M[k, j] = mutual(coils[k], coils[j]; xorder=order, yorder=order)
            (j != k) && (M[j, k] = M[k, j])
        end
    end

    Mstar = deepcopy(M)
    for j in eachindex(coils)
        for k in eachindex(coils)[j:end]
            Mstar[k, j] -= b[k] * b[j] / K
            (j != k) && (Mstar[j, k] = Mstar[k, j])
        end
    end

    # reuse b vector for resistances
    for j in eachindex(coils)
        b[j] = coils[j].resistance
    end
    R = Diagonal(b)

    A0 = .- Mstar \ R

    D, V = eigen(A0)
    dmax = argmax(real.(D))
    γ = real(D[dmax])
    passive = [current(C) == 0.0 for C in coils]
    @views v = V[passive, dmax]
    @views τ = dot(v, M[passive, passive], v) / dot(v, R[passive, passive], v)
    return γ, τ, γ * τ
end


