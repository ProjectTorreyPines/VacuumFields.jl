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

function Image(EQ::MXHEquilibrium.AbstractEquilibrium)
    Sb = MXHEquilibrium.plasma_boundary(EQ; precision=0.0)
    if Sb === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        Sb = MXHEquilibrium.plasma_boundary(EQ)
    end
    return Image(EQ, Sb.r, Sb.z)
end

function Image(EQ::MXHEquilibrium.AbstractEquilibrium, Rb::AbstractVector{T}, Zb::AbstractVector{T}) where {T<:Real}
    Lb = cumlength(Rb, Zb)
    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol
    cocos = MXHEquilibrium.cocos(EQ)
    fac = cocos.sigma_rhotp * cocos.sigma_Bp * (2π)^cocos.exp_Bp
    dψdn_R = similar(Lb)
    dψdn_R .= fac .* MXHEquilibrium.poloidal_Bfield.(Ref(EQ), Rb, Zb)
    return Image(Rb, Zb, Lb, dψdn_R)
end

function Image(dd::IMAS.dd)
    return Image(dd.equilibrium.time_slice[])
end

function Image(eqt::IMAS.equilibrium__time_slice)
    Rb = eqt.boundary.outline.r
    Zb = eqt.boundary.outline.z

    Lb = cumlength(Rb, Zb)
    # dPsi/dn = σ_RφZ * σ_ρθφ * Bp_fac * R * Bpol

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt2d)
    # poloidal magnetic field (with sign)
    Br, Bz = IMAS.Br_Bz(PSI_interpolant, Rb, Zb)
    dψdn_R = (sqrt.(Br .^ 2.0 .+ Bz .^ 2.0) .*
              sign.((Zb .- eqt.global_quantities.magnetic_axis.z) .* Br .- (Rb .- eqt.global_quantities.magnetic_axis.r) .* Bz)
    )
    dψdn_R .*= 2π # σ_RφZ * σ_ρθφ * Bp_fac in COCOS 11
    return Image(Rb, Zb, Lb, dψdn_R)
end