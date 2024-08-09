# ******
# Utils
# ******
function cost_λ_regularize(
    λ_exp::Float64,
    coils::Vector{<:AbstractCoil},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

    _, cost = find_coil_currents!(coils, EQ, image, flux_cps, saddle_cps; ψbound, fixed_coils, λ_regularize=10^λ_exp, Sb)

    return cost^2
end

"""
    optimal_λ_regularize(
        coils::Vector{<:AbstractCoil},
        EQ::MXHEquilibrium.AbstractEquilibrium,
        image::Image,
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
        ψbound::Real=0.0,
        fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        min_exp::Integer=-20, max_exp::Integer=-10,
        Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

Find the optimizal `λ_regularize` to be used within `find_coil_currents!` that minimizes the cost
"""
function optimal_λ_regularize(
    coils::Vector{<:AbstractCoil},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    min_exp::Integer=-20, max_exp::Integer=-10,
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

    λ_range_exp = min_exp:0.5:max_exp
    cost_λ = λ -> cost_λ_regularize(λ, coils, EQ, image, flux_cps, saddle_cps; ψbound, fixed_coils, Sb)
    return 10^λ_range_exp[argmin(cost_λ(λ) for λ in λ_range_exp)]
end

# ***************************************************
# Transform fixed-boundary ψ to free-boundary ψ
# ***************************************************
"""
    encircling_fixed2free(
        EQfixed::MXHEquilibrium.AbstractEquilibrium,
        n_coils::Integer;
        Rx::AbstractVector{Float64}=Float64[],
        Zx::AbstractVector{Float64}=Float64[],
        fraction_inside::Float64=1.0 - 1E-4,
        λ_regularize::Float64=-1.0,
        Rgrid::AbstractVector{Float64}=EQfixed.r,
        Zgrid::AbstractVector{Float64}=EQfixed.z)

Distribute n point coils around fixed boundary plasma to get a free boundary ψ map
"""
function encircling_fixed2free(
    shot::TEQUILA.Shot,
    n_coils::Integer;
    n_grid=129,
    scale::Float64=1.2,
    kwargs...)

    R0 = shot.surfaces[1, end]
    Z0 = shot.surfaces[2, end]
    ϵ = shot.surfaces[3, end]
    κ = shot.surfaces[4, end]
    a = R0 * ϵ * scale
    b = κ * a

    Rgrid = range(max(R0 - a, 0.0), R0 + a, n_grid)
    Zgrid = range(Z0 - b, Z0 + b, n_grid)
    return Rgrid, Zgrid, encircling_fixed2free(shot, n_coils, Rgrid, Zgrid; kwargs...)
end

function encircling_fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    n_coils::Integer,
    Rgrid::AbstractVector{Float64}=EQfixed.r,
    Zgrid::AbstractVector{Float64}=EQfixed.z;
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    fraction_inside::Float64=0.999,
    λ_regularize::Float64=0.0,
    ψbound::Float64=0.0)

    image = Image(EQfixed)
    coils = encircling_coils(EQfixed, n_coils)
    if λ_regularize < 0.0
        λ_regularize = optimal_λ_regularize(coils, EQfixed, image, flux_cps, saddle_cps; ψbound)
    end
    append!(flux_cps, boundary_control_points(EQfixed, fraction_inside, ψbound))

    find_coil_currents!(coils, EQfixed, image, flux_cps, saddle_cps; ψbound, λ_regularize)

    return fixed2free(EQfixed, coils, Rgrid, Zgrid)
end

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

plasma_boundary_psi_w_fallback(shot::TEQUILA.Shot, args...) = MXHEquilibrium.Boundary(MXHEquilibrium.MXH(shot, 1.0)()...), 0.0

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
        ψbound::Real=0.0,
        fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Real=0.0) where {T<:Real}

Convert the flux of a fixed-boundary equilibrium `EQfixed` to a free-boundary representation on an `(R,Z)` grid,
    subtracting out the flux from image currents `image` and adding in the flux from `coils` with currents
    that best satisfy the flux and saddle control points and `fixed_coils` with predefined coil currents
"""
function fixed2free(
    EQfixed::MXHEquilibrium.AbstractEquilibrium,
    image::Image,
    coils::AbstractVector{<:AbstractCoil},
    R::AbstractVector{T},
    Z::AbstractVector{T};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Real=0.0) where {T<:Real}

    Sb, ψb = plasma_boundary_psi_w_fallback(EQfixed)
    ψ_f2f = T[MXHEquilibrium.in_boundary(Sb, (r, z)) ? EQfixed(r, z) - ψb + ψbound : ψbound for z in Z, r in R]

    if (!isempty(flux_cps) || !isempty(saddle_cps))
        if λ_regularize < 0.0
            λ_regularize = optimal_λ_regularize(coils, EQfixed, image, flux_cps, saddle_cps; ψbound, fixed_coils, Sb)
        end
        find_coil_currents!(coils, EQfixed, image, flux_cps, saddle_cps; ψbound, fixed_coils, λ_regularize, Sb)
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
    coils::Vector{<:AbstractCoil},
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
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

    @series begin
        seriestype --> :contour
        levels --> levels
        color --> :black
        linewidth --> 2
        R, Z, free
    end

end