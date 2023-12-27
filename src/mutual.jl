# mutual inductance between coils

function mutual(C1::AbstractCoil, C2::PointCoil;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    fac = μ₀ * Bp_fac * C1.turns * C2.turns
    return fac * Green(C1, C2.R, C2.Z)
end

function mutual(C1::AbstractCoil, C2::DistributedCoil;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    fac = μ₀ * Bp_fac * C1.turns * C2.turns
    return fac * sum(Green(C1, C2.R[k], C2.Z[k]) for k in eachindex(C2.R)) / length(C2.R)
end

function mutual(C1::AbstractCoil, C2::ParallelogramCoil;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp,
                xorder::Int=3, yorder::Int=3)
    fac = μ₀ * Bp_fac * C1.turns * C2.turns
    f = (r, z) -> Green(C1, r, z; xorder = xorder+1, yorder = yorder+1)
    return fac * integrate(f, C2; xorder, yorder) / area(C2)
end


# Plasma-coil mutuals
function mutual(EQ::MXHEquilibrium.AbstractEquilibrium, C::AbstractCoil, δZ::Real=0.0; kwargs...)
    return mutual(EQ, Image(EQ), C, δZ; kwargs...)
end

function dM_dZ(EQ::MXHEquilibrium.AbstractEquilibrium, C::AbstractCoil; kwargs...)
    return dM_dZ(EQ, Image(EQ), C; kwargs...)
end

# Shifting plasma up by δZ is the same as shifting the coil down by δZ
_pfunc(Pfunc, image::Image, C::PointCoil, δZ; kwargs...) = Pfunc(image, C.r, C.Z - δZ)

function _pfunc(Pfunc, image::Image, C::DistributedCoil, δZ; kwargs...)
    return sum(Pfunc(image, C.R[k], C.Z[k] - δZ) for k in eachindex(C.R)) / length(C.R)
end

function _pfunc(Pfunc, image::Image, C::ParallelogramCoil, δZ; xorder::Int=3, yorder::Int=3, kwargs...)
    f = (r, z) -> Pfunc(image, r, z - δZ)
    return integrate(f, C; xorder, yorder) / area(C)
end

function mutual(EQ::MXHEquilibrium.AbstractEquilibrium, image::Image, C::AbstractCoil, δZ::Real=0.0)
    Ψ = _pfunc(ψ, image, C, δZ)
    Ip = MXHEquilibrium.plasma_current(EQ)
    return C.turns * Ψ / Ip
end

function dM_dZ(EQ::MXHEquilibrium.AbstractEquilibrium, image::Image, C::AbstractCoil, δZ::Real=0.0)
    # dψ/d(δZ) = -dψ_dZ
    Ψ = -_pfunc(dψ_dZ, image, C, δZ)
    Ip = MXHEquilibrium.plasma_current(EQ)
    return C.turns * Ψ / Ip
end

function d2M_dZ2(EQ::MXHEquilibrium.AbstractEquilibrium, image::Image, C::AbstractCoil, δZ::Real=0.0)
    # dψ/d(δZ) = -dψ_dZ
    Ψ = _pfunc(d2ψ_dZ2, image, C, δZ)
    Ip = MXHEquilibrium.plasma_current(EQ)
    return C.turns * Ψ / Ip
end