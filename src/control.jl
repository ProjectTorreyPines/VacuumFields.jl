abstract type AbstractControlPoint{T<:Real} end

mutable struct FluxControlPoint{T<:Real} <: AbstractControlPoint{T}
    R::T
    Z::T
    target::T
    weight::T
end

FluxControlPoint(R, Z, target) = FluxControlPoint(R, Z, target, one(R))
FluxControlPoint(R, Z, target, weight) = FluxControlPoint(promote(R, Z, target, weight)...)

mutable struct SaddleControlPoint{T<:Real} <: AbstractControlPoint{T}
    R::T
    Z::T
    weight::T
end

SaddleControlPoint(R, Z) = SaddleControlPoint(R, Z, one(R))
SaddleControlPoint(R, Z, weight) = SaddleControlPoint(promote(R, Z, weight)...)

@recipe function plot_SaddleControlPoint(cp::SaddleControlPoint)
    @series begin
        marker := :star
        cp, nothing
    end
end

@recipe function plot_FluxControlPoint(cp::FluxControlPoint)
    @series begin
        marker := :circle
        cp, nothing
    end
end

@recipe function plot_control_point(cp::AbstractControlPoint, dispatch::Nothing=nothing)
    @series begin
        seriestype := :scatter
        label --> ""
        [cp.R], [cp.Z]
    end
end

@recipe function plot_control_points(cps::AbstractVector{<:AbstractControlPoint})
    for (k, cp) in enumerate(cps)
        @series begin
            primary --> k == 1
            cp
        end
    end
end

function reg_solve(A, b, λ)
    return (A' * A + λ * I) \ A' * b
end

function FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtarget::Real=0.0)
    return [FluxControlPoint(Rs[k], Zs[k], ψtarget) for k in eachindex(Rs)]
end
function FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtarget::AbstractVector{<:Real})
    return [FluxControlPoint(Rs[k], Zs[k], ψtarget[k]) for k in eachindex(Rs)]
end

function SaddleControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real})
    return [SaddleControlPoint(Rs[k], Zs[k]) for k in eachindex(Rs)]
end

function boundary_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
    ψ0, ψb = MXHEquilibrium.psi_limits(EQfixed)
    _, ψb = plasma_boundary_psi_w_fallback(EQfixed)

    Sp = MXHEquilibrium.flux_surface(EQfixed, fraction_inside * (ψb - ψ0) + ψ0; n_interp=Npts)
    ψtarget = fraction_inside * (ψb - ψ0) + ψ0 - ψb + ψbound
    return [FluxControlPoint(Sp.r[k], Sp.z[k], ψtarget) for k in eachindex(Sp.r)]
end

function boundary_control_points(shot::TEQUILA.Shot, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
    bnd = MillerExtendedHarmonic.MXH(shot, fraction_inside)
    θs = LinRange(0, 2π, Npts + 1)
    ψtarget = ψbound + TEQUILA.psi_ρθ(shot, fraction_inside, 0.0)
    return [FluxControlPoint(bnd(θ)..., ψtarget) for θ in θs[1:end-1]]
end

function find_coil_currents!(
    coils::Vector{<:AbstractCoil},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0)

    return find_coil_currents!(coils, EQ, Image(EQ), flux_cps, saddle_cps; ψbound, fixed_coils, λ_regularize)
end

function find_coil_currents!(
    coils::Vector{<:AbstractCoil},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0)

    # First reset current in coils to unity
    for coil in coils
        coil.current = 1.0
    end

    N = length(flux_cps) + 2 * length(saddle_cps)
    A = zeros(N, length(coils))
    b = zeros(N)

    init_b!(b, EQ, image, flux_cps, saddle_cps; ψbound)
    cocos = MXHEquilibrium.cocos(EQ)
    populate_Ab!(A, b, coils, flux_cps, saddle_cps; fixed_coils, cocos)

    # Least-squares solve for coil currents
    if λ_regularize > 0
        # Least-squares with regularization
        # https://www.youtube.com/watch?v=9BckeGN0sF0
        Ic0 = reg_solve(A, b, λ_regularize / length(coils)^2)
    else
        Ic0 = A \ b
    end

    # update values of coils current
    for (k, coil) in enumerate(coils)
        coil.current = Ic0[k]
    end

    cost = norm(A * Ic0 .- b) / norm(b)

    return Ic0, cost
end

function find_coil_currents!(
    coils::Vector{<:AbstractCoil},
    EQ::Nothing,  # VACUUM case
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11))

    return find_coil_currents!(coils, flux_cps, saddle_cps; fixed_coils, λ_regularize, cocos)
end

function find_coil_currents!(
    coils::Vector{<:AbstractCoil},
    EQ::Nothing, # VACUUM case
    image::Nothing,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0,
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11))

    return find_coil_currents!(coils, flux_cps, saddle_cps; fixed_coils, λ_regularize, cocos)
end

function find_coil_currents!(
    coils::Vector{<:AbstractCoil},
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11))

    # First reset current in coils to unity
    for coil in coils
        coil.current = 1.0
    end

    N = length(flux_cps) + 2 * length(saddle_cps)
    A = zeros(N, length(coils))
    b = zeros(N)

    populate_Ab!(A, b, coils, flux_cps, saddle_cps; fixed_coils, cocos)

    # Least-squares solve for coil currents
    if λ_regularize > 0
        # Least-squares with regularization
        # https://www.youtube.com/watch?v=9BckeGN0sF0
        Ic0 = reg_solve(A, b, λ_regularize / length(coils)^2)
    else
        Ic0 = A \ b
    end

    # update values of coils current
    for (k, coil) in enumerate(coils)
        coil.current = Ic0[k]
    end

    cost = norm(A * Ic0 .- b) / norm(b)

    return Ic0, cost
end

function init_b!(
    b::AbstractVector{<:Real},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    ψbound::Real=0.0)

    Nflux = length(flux_cps)
    _, ψb = MXHEquilibrium.psi_limits(EQ)

    @threads for i in eachindex(flux_cps)
        cp = flux_cps[i]
        r = cp.R
        z = cp.Z

        # RHS

        # remove plasma current contribution (EQ - image)
        b[i] -= ψbound - ψb + EQ(r, z)
        b[i] += ψ(image, r, z)
    end

    @threads for i in eachindex(saddle_cps)
        cp = saddle_cps[i]
        r = cp.R
        z = cp.Z

        ir = Nflux + 2i - 1
        iz = Nflux + 2i

        # target is zero and assume outside boundary, so no EQ contribution
        b[ir] += dψ_dR(image, r, z)
        b[iz] += dψ_dZ(image, r, z)
    end

    return b
end

function populate_Ab!(A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real},
    coils::Vector{<:AbstractCoil},
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[];
    fixed_coils::Vector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    cocos=MXHEquilibrium.cocos(11))

    Nflux = length(flux_cps)

    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp

    @threads for i in eachindex(flux_cps)
        cp = flux_cps[i]
        r = cp.R
        z = cp.Z

        # RHS

        # target
        b[i] += cp.target

        # remove fixed coil contribution
        if !isempty(fixed_coils)
            b[i] -= sum(ψ(fixed_coil, r, z; Bp_fac) for fixed_coil in fixed_coils)
        end

        # Build matrix relating coil Green's functions to boundary points
        A[i, :] .= ψ.(coils, r, z; Bp_fac)

        # weighting
        w = sqrt(cp.weight)
        b[i] *= w
        A[i, :] .*= w
    end

    @threads for i in eachindex(saddle_cps)
        cp = saddle_cps[i]
        r = cp.R
        z = cp.Z

        ir = Nflux + 2i - 1
        iz = Nflux + 2i

        # remove fixed coil contribution
        if !isempty(fixed_coils)
            b[ir] -= sum(dψ_dR(fixed_coil, r, z; Bp_fac) for fixed_coil in fixed_coils)
            b[iz] -= sum(dψ_dZ(fixed_coil, r, z; Bp_fac) for fixed_coil in fixed_coils)
        end

        # Build matrix relating coil Green's functions to boundary points
        A[ir, :] .= dψ_dR.(coils, r, z; Bp_fac)
        A[iz, :] .= dψ_dZ.(coils, r, z; Bp_fac)

        # weighting
        w = sqrt(cp.weight)
        b[ir:iz] .*= w
        A[ir:iz, :] .*= w
    end
end