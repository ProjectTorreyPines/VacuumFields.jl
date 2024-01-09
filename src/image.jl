struct Image{V<:AbstractVector{<:Real}}
    Rb::V
    Zb::V
    Lb::V
    dψdn_R::V
    _Vb::Vector{V}
end

function Image(Rb, Zb, Lb, dψdn_R)
    Vb = [zero(Lb) for _ in 1:Threads.nthreads()]
    return Image(Rb, Zb, Lb, dψdn_R, Vb)
end

function cumlength(R::T, Z::T) where {T<:AbstractVector{Float64}}
    # Length along boundary
    L = Array{Float64,1}(undef, length(R))
    @inbounds L[1] = 0
    for i in eachindex(R)[2:end]
        @inbounds L[i] = L[i-1] + sqrt((R[i] - R[i-1])^2 + (Z[i] - Z[i-1])^2)
    end
    return L
end

function Image(shot::TEQUILA.Shot)
    bnd = @views MillerExtendedHarmonic.MXH(shot.surfaces[:, end])
    return Image(shot, bnd)
end
function Image(shot::TEQUILA.Shot, bnd::MillerExtendedHarmonic.MXH; Nb::Integer=100 * length(bnd.c))
    image = Image(shot, bnd(Nb; adaptive=false)...)
    return image
end

function Image(EQ::MXHEquilibrium.AbstractEquilibrium)
    _, Sb = MXHEquilibrium.plasma_boundary_psi(EQ; precision=0.0)
    if Sb === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        _, Sb = MXHEquilibrium.plasma_boundary_psi(EQ)
    end
    return Image(EQ, Sb.r, Sb.z)
end

function Image(EQ::MXHEquilibrium.AbstractEquilibrium, Rb::AbstractVector{T}, Zb::AbstractVector{T}) where {T<:Real}
    Lb = cumlength(Rb, Zb)
    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    cocos = MXHEquilibrium.cocos(EQ)
    fac = cocos.sigma_rhotp * cocos.sigma_Bp * (2π)^cocos.exp_Bp
    dψdn_R = zero(Lb)
    dψdn_R .= fac .* MXHEquilibrium.poloidal_Bfield.(Ref(EQ), Rb, Zb)
    return Image(Rb, Zb, Lb, dψdn_R)
end