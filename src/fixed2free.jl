# ***************************************************
# Transform fixed-boundary ψ to free-boundary ψ given currents in coils
# ***************************************************
function plasma_boundary_psi_w_fallback(EQ::MXHEquilibrium.AbstractEquilibrium)
    Sb, ψb = MXHEquilibrium.plasma_boundary_psi(EQ; precision=0.0)
    if Sb === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        Sb, ψb = MXHEquilibrium.plasma_boundary(EQ)
    end
    return Sb, ψb
end

"""
    fixed2free(
        EQfixed::MXHEquilibrium.AbstractEquilibrium,
        coils::AbstractVector{<:AbstractCoil},
        R::AbstractVector{T},
        Z::AbstractVector{T};
        ψbound::Real=0.0,
        kwargs...) where {T<:Real}

Convert the flux of a fixed-boundary equilibrium `EQfixed` to a free-boundary representation on an `(R,Z)` grid,
using the flux from `coils` with currents satisfying given control points
"""
function fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    coils::AbstractVector{<:AbstractCoil},
    R::AbstractVector{T},
    Z::AbstractVector{T};
    ψbound::Real=0.0,
    kwargs...) where {T<:Real}
    image = Image(EQfixed)
    return fixed2free(EQfixed, image, coils, R, Z; ψbound, kwargs...)
end

"""
    fixed2free(
        EQfixed::MXHEquilibrium.AbstractEquilibrium,
        image::Image,
        coils::AbstractVector{<:AbstractCoil},
        R::AbstractVector{T},
        Z::AbstractVector{T};
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Real=0.0) where {T<:Real}

Convert the flux of a fixed-boundary equilibrium `EQfixed` to a free-boundary representation on an `(R,Z)` grid,
subtracting out the flux from image currents `image` and adding in the flux from `coils` with currents
that best satisfy the control points and `fixed_coils` with predefined coil currents
"""
function fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    image::Image,
    coils::AbstractVector{<:AbstractCoil},
    R::AbstractVector{T},
    Z::AbstractVector{T};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Real=0.0) where {T<:Real}

    Sb, ψb = plasma_boundary_psi_w_fallback(EQfixed)
    ψ_f2f = T[MXHEquilibrium.in_boundary(Sb, (r, z)) ? EQfixed(r, z) - ψb + ψbound : ψbound for z in Z, r in R]

    if (!isempty(flux_cps) || !isempty(saddle_cps) || !isempty(iso_cps))
        find_coil_currents!(coils, EQfixed, image; flux_cps, saddle_cps, iso_cps, ψbound, fixed_coils, λ_regularize, Sb)
    end

    # ψ from image and coil currents
    cocos = MXHEquilibrium.cocos(EQfixed)
    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp
    Threads.@threads for i in eachindex(R)
        @inbounds r = R[i]
        for j in eachindex(Z)
            @inbounds z = Z[j]
            # subtract image ψ
            @inbounds ψ_f2f[j, i] -= ψ(image, r, z)
            # add coil ψ
            @inbounds ψ_f2f[j, i] += sum(ψ(coil, r, z; Bp_fac) for coil in coils)
        end
    end

    # sometimes ψ_f2f goes to Inf
    ψ_f2f[abs.(ψ_f2f).==Inf] .= 0.0

    return ψ_f2f
end

@recipe function plot_free(
    R::AbstractVecOrMat{<:Real},
    Z::AbstractVecOrMat{<:Real},
    free::AbstractVecOrMat{<:Real},
    coils::AbstractVector{<:AbstractCoil};
    flux_cps=FluxControlPoint{Float64}[],
    saddle_cps=SaddleControlPoint{Float64}[],
    iso_cps=IsoControlPoint{Float64}[],
    clims=(-1, 1) .* maximum(abs, free))

    aspect_ratio --> :equal
    xlims --> extrema(R)
    ylims --> extrema(Z)

    levels = Float64[]

    @series begin
        seriestype --> :heatmap
        c --> :diverging
        clims --> clims
        R, Z, free
    end

    for (k, coil) in enumerate(coils)
        @series begin
            label --> ((k == 1) ? "Coil" : :none)
            coil
        end
    end

    for (k, cp) in enumerate(flux_cps)
        !(cp.target in levels) && push!(levels, cp.target)
        @series begin
            seriestype --> :scatter
            color --> :cyan
            label --> ((k == 1) ? "Flux CP" : :none)
            [cp.R], [cp.Z]
        end
    end

    for (k, cp) in enumerate(saddle_cps)
        @series begin
            seriestype --> :scatter
            color --> :magenta
            label --> ((k == 1) ? "Saddle CP" : :none)
            [cp.R], [cp.Z]
        end
    end

    for (k, cp) in enumerate(iso_cps)
        @series begin
            seriestype --> :scatter
            color --> :yellow
            label --> ((k == 1) ? "Iso CP" : :none)
            [cp.R1, cp.R2], [cp.Z1, cp.Z2]
        end
    end

    @series begin
        seriestype --> :contour
        levels --> levels
        color --> :black
        linewidth --> 2
        R, Z, free
    end

end