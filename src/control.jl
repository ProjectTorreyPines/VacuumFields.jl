"""
    AbstractControlPoint{T<:Real}

Abstract control point - for coil current least-squres fitting
"""
abstract type AbstractControlPoint{T<:Real} end

mutable struct FluxControlPoint{T<:Real} <: AbstractControlPoint{T}
    R::T
    Z::T
    target::T
    weight::T
end

"""
    FluxControlPoint(R::Real, Z::Real, target::Real, weight::Real=1.0)

Return a control point for a `target` flux value at point `(R, Z)`, with an optional `weight`
"""
FluxControlPoint(R::Real, Z::Real, target::Real, weight::Real=1.0) = FluxControlPoint(promote(R, Z, target, weight)...)

mutable struct SaddleControlPoint{T<:Real} <: AbstractControlPoint{T}
    R::T
    Z::T
    weight::T
end

"""
    SaddleControlPoint(R::Real, Z::Real,weight::Real=1.0)

Return a control point for a saddle (i.e., dψ/dR = dψ/dZ = 0.0) at point `(R, Z)`, with an optional `weight`
"""
SaddleControlPoint(R::Real, Z::Real, weight::Real=1.0) = SaddleControlPoint(promote(R, Z, weight)...)


mutable struct IsoControlPoint{T<:Real} <: AbstractControlPoint{T}
    R1::T
    Z1::T
    R2::T
    Z2::T
    weight::T
end

"""
    IsoControlPoint(R::Real, Z::Real,weight::Real=1.0)

Return a control point for equal flux between points `(R1, Z1)` and `(R2, Z2)`, with an optional `weight`
"""
IsoControlPoint(R1::Real, Z1::Real, R2::Real, Z2::Real, weight::Real=1.0) = IsoControlPoint(promote(R1, Z1, R2, Z2, weight)...)


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

@recipe function plot_IsoControlPoint(cp::IsoControlPoint)
    @series begin
        marker := :diamond
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

"""
    FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtarget::Real)

Return a Vector of FluxControlPoint at each `Rs` and `Zs` point, each with the same `ψtarget` flux
"""
function FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtarget::Real)
    return [FluxControlPoint(Rs[k], Zs[k], ψtarget) for k in eachindex(Rs)]
end

"""
    FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtargets::AbstractVector{<:Real})

Return a Vector of FluxControlPoint at each `Rs` and `Zs` point with `ψtargets` flux
"""
function FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtargets::AbstractVector{<:Real})
    return [FluxControlPoint(Rs[k], Zs[k], ψtargets[k]) for k in eachindex(Rs)]
end

"""
    SaddleControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real})

Return a Vector of SaddleControlPoint at each `Rs` and `Zs` point
"""
function SaddleControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real})
    return [SaddleControlPoint(Rs[k], Zs[k]) for k in eachindex(Rs)]
end

"""
    IsoControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real})

Return a Vector of IsoControlPoints between each pair of `Rs` and `Zs` points
"""
function IsoControlPoints(Rs::AbstractVector{T}, Zs::AbstractVector{T}) where {T <: Real}
    N = length(Rs)
    @assert length(Zs) == N
    is_closed = (Rs[1] ≈ Rs[end]) && (Zs[1] ≈ Zs[end])
    Nmax = is_closed ? N-1 : N
    iso_cps = [IsoControlPoint(Rs[k], Zs[k], Rs[1], Zs[1]) for k in eachindex(Rs)[2:Nmax]]
    return iso_cps
end

"""
    boundary_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
Return a Vector of FluxControlPoint, each with target `ψbound`, at `Npts` equally distributed `fraction_inside` percent inside the the boundary of `EQfixed`
"""
function boundary_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
    ψ0, ψb = MXHEquilibrium.psi_limits(EQfixed)
    Sp = MXHEquilibrium.flux_surface(EQfixed, fraction_inside * (ψb - ψ0) + ψ0; n_interp=Npts)
    ψtarget = fraction_inside * (ψb - ψ0) + ψ0 - ψb + ψbound
    return [FluxControlPoint(Sp.r[k], Sp.z[k], ψtarget, 1.0 / Npts) for k in 1:length(Sp.r)-1]
end

"""
    boundary_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
Return a Vector of FluxControlPoint, each with target `ψbound`, at `Npts` equally distributed `fraction_inside` percent inside the the boundary of `shot`
"""
function boundary_control_points(shot::TEQUILA.Shot, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
    bnd = MillerExtendedHarmonic.MXH(shot, fraction_inside)
    θs = LinRange(0, 2π, Npts + 1)
    ψtarget = ψbound + TEQUILA.psi_ρθ(shot, fraction_inside, 0.0)
    return [FluxControlPoint(bnd(θ)..., ψtarget, 1.0 / Npts) for θ in θs[1:end-1]]
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        EQ::MXHEquilibrium.AbstractEquilibrium;
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

Find the currents for `coils` that best match (least-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Assumes flux from  plasma current given by equilibrium `EQ` with a `ψbound` flux at the boundary `Sb`

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    EQ::MXHEquilibrium.AbstractEquilibrium;
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

    return find_coil_currents!(coils, EQ, Image(EQ); flux_cps, saddle_cps, iso_cps, ψbound, fixed_coils, λ_regularize, Sb)
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        EQ::MXHEquilibrium.AbstractEquilibrium,
        image::Image;
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Assumes flux from  plasma current given by equilibrium `EQ` with image currents `image` and a `ψbound` flux at the boundary `Sb`

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image;
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

    # First reset current in coils to unity
    for coil in coils
        set_current!(coil, 1.0)
    end

    N = length(flux_cps) + 2 * length(saddle_cps) + length(iso_cps)
    A = zeros(N, length(coils))
    b = zeros(N)

    init_b!(b, EQ, image; flux_cps, saddle_cps, iso_cps, ψbound, Sb)
    cocos = MXHEquilibrium.cocos(EQ)
    populate_Ab!(A, b, coils; flux_cps, saddle_cps, iso_cps, fixed_coils, cocos)

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
        set_current!(coil, Ic0[k])
    end

    cost = norm(A * Ic0 .- b) / norm(b)

    return Ic0, cost
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        EQ::MXHEquilibrium.AbstractEquilibrium,
        image::Image;
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Assumes flux from  plasma current given by equilibrium `EQ` with image currents `image` and a `ψbound` flux at the boundary `Sb`

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    ψpl::Interpolations.AbstractInterpolation,
    dψpl_dR::Union{Function, Interpolations.AbstractInterpolation} = (r, z) -> Interpolations.gradient(ψpl, r, z)[1],
    dψpl_dZ::Union{Function, Interpolations.AbstractInterpolation} = (r, z) -> Interpolations.gradient(ψpl, r, z)[2];
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11),
    Sb=nothing)

    # First reset current in coils to unity
    for coil in coils
        set_current!(coil, 1.0)
    end

    N = length(flux_cps) + 2 * length(saddle_cps) + length(iso_cps)
    A = zeros(N, length(coils))
    b = zeros(N)

    init_b!(b, ψpl, dψpl_dR, dψpl_dZ; flux_cps, saddle_cps, iso_cps)
    populate_Ab!(A, b, coils; flux_cps, saddle_cps, iso_cps, fixed_coils, cocos)

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
        set_current!(coil, Ic0[k])
    end

    cost = norm(A * Ic0 .- b) / norm(b)

    return Ic0, cost
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        EQ::Nothing;
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        cocos=MXHEquilibrium.cocos(11),
        Sb=nothing)

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Vacuume case: assumes no equilibrium plasma current

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    EQ::Nothing;  # VACUUM case
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11),
    Sb=nothing)

    return find_coil_currents!(coils; flux_cps, saddle_cps, iso_cps, fixed_coils, λ_regularize, cocos)
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        EQ::Nothing, # VACUUM case
        image::Nothing;
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        cocos=MXHEquilibrium.cocos(11),
        Sb=nothing)

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Vacuume case: assumes no equilibrium plasma current

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    EQ::Nothing, # VACUUM case
    image::Nothing;
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11),
    Sb=nothing)

    return find_coil_currents!(coils; flux_cps, saddle_cps, iso_cps, fixed_coils, λ_regularize, cocos)
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil};
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        cocos=MXHEquilibrium.cocos(11))

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Vacuume case: assumes no equilibrium plasma current

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11))

    # First reset current in coils to unity
    for coil in coils
        set_current!(coil, 1.0)
    end

    N = length(flux_cps) + 2 * length(saddle_cps) + length(iso_cps)
    A = zeros(N, length(coils))
    b = zeros(N)

    populate_Ab!(A, b, coils; flux_cps, saddle_cps, iso_cps, fixed_coils, cocos)

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
        set_current!(coil, Ic0[k])
    end

    cost = norm(A * Ic0 .- b) / norm(b)

    return Ic0, cost
end

function init_b!(
    b::AbstractVector{T},
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image;
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1]) where {T <: Real}

    _, ψb = MXHEquilibrium.psi_limits(EQ)

    ψpl  = (r, z) -> (MXHEquilibrium.in_boundary(Sb, (r, z)) ? EQ(r, z) : ψb) +  ψbound - ψb - ψ(image, r, z)
    dψpl_dR = (r, z) -> -dψ_dR(image, r, z)
    dψpl_dZ = (r, z) -> -dψ_dZ(image, r, z)

    init_b!(b, ψpl, dψpl_dR, dψpl_dZ; flux_cps, saddle_cps, iso_cps)

end

function init_b!(
    b::AbstractVector{T},
    ψpl::Union{Function, Interpolations.AbstractInterpolation},
    dψpl_dR::Union{Function, Interpolations.AbstractInterpolation},
    dψpl_dZ::Union{Function, Interpolations.AbstractInterpolation};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[]) where {T <: Real}


    b .= 0.0

    Nflux = length(flux_cps)
    Nsaddle = length(saddle_cps)

    @threads for i in eachindex(flux_cps)
        cp = flux_cps[i]
        r = cp.R
        z = cp.Z

        # subtract plasma current contribution
        b[i] -= ψpl(r, z)
    end

    @threads for i in eachindex(saddle_cps)
        cp = saddle_cps[i]
        r = cp.R
        z = cp.Z

        ir = Nflux + 2i - 1
        iz = Nflux + 2i

        # subtract plasma contribution
        b[ir] -= dψpl_dR(r, z)
        b[iz] -= dψpl_dZ(r, z)
    end

    r1_old, z1_old, r2_old, z2_old = -one(T), zero(T), -one(T), zero(T)
    ψ1_old, ψ2_old = zero(T), zero(T)
    @threads for i in eachindex(iso_cps)
        cp = iso_cps[i]
        r1, z1 = cp.R1, cp.Z1
        r2, z2 = cp.R2, cp.Z2

        # check if (r1, z1) was used in last iteration
        if (r1 == r2_old) && (z1 == z2_old)
            ψ1 = ψ2_old
        elseif (r1 == r1_old) && (z1 == z1_old)
            ψ1 = ψ1_old
        else
            ψ1 = ψpl(r1, z1)
        end

        # check if (r2, z2) was used in last iteration
        if (r2 == r1_old) && (z2 == z1_old)
            ψ2 = ψ1_old
        elseif (r2 == r2_old) && (z2 == z2_old)
            ψ2 = ψ2_old
        else
            ψ2 = ψpl(r2, z2)
        end

        k = Nflux + 2Nsaddle + i
        b[k] -= ψ1 - ψ2

        # store for next iteration
        r1_old, z1_old, r2_old, z2_old = r1, z1, r2, z2
        ψ1_old, ψ2_old = ψ1, ψ2
    end

    return b
end

function populate_Ab!(A::AbstractMatrix{T}, b::AbstractVector{T},
    coils::AbstractVector{<:AbstractCoil};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    cocos=MXHEquilibrium.cocos(11)) where {T <: Real}

    Nflux = length(flux_cps)
    Nsaddle = length(saddle_cps)

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

    if !isempty(iso_cps)
        r1_old, z1_old, r2_old, z2_old = -one(T), zero(T), -one(T), zero(T)
        ψf1_old, ψf2_old = zero(T), zero(T)
        Ncoil = length(coils)
        ψc1_old, ψc2_old = Vector{T}(undef, Ncoil), Vector{T}(undef, Ncoil)
        ψc1, ψc2 = Vector{T}(undef, Ncoil), Vector{T}(undef, Ncoil)
        for i in eachindex(iso_cps)
            cp = iso_cps[i]
            r1, z1 = cp.R1, cp.Z1
            r2, z2 = cp.R2, cp.Z2

            k = Nflux + 2Nsaddle + i

            # remove fixed coil contribution
            if !isempty(fixed_coils)
                if (r1 == r2_old) && (z1 == z2_old)
                    ψf1 = ψf2_old
                elseif (r1 == r1_old) && (z1 == z1_old)
                    ψf1 = ψf1_old
                else
                    ψf1 = sum(ψ(fixed_coil, r1, z1; Bp_fac) for fixed_coil in fixed_coils)
                end

                if (r2 == r1_old) && (z2 == z1_old)
                    ψf2 = ψf1_old
                elseif (r2 == r2_old) && (z2 == z2_old)
                    ψf2 = ψf2_old
                else
                    ψf2 = sum(ψ(fixed_coil, r2, z2; Bp_fac) for fixed_coil in fixed_coils)
                end

                b[k] -= ψf1 - ψf2

                ψf1_old, ψf2_old = ψf1, ψf2
            end

            # Build matrix relating coil Green's functions to boundary points

            if (r1 == r2_old) && (z1 == z2_old)
                ψc1 .= ψc2_old
            elseif (r1 == r1_old) && (z1 == z1_old)
                ψc1 .= ψc1_old
            else
                ψc1 .= ψ.(coils, r1, z1; Bp_fac)
            end

            if (r2 == r1_old) && (z2 == z1_old)
                ψc2 .= ψc1_old
            elseif (r2 == r2_old) && (z2 == z2_old)
                ψc2 .= ψc2_old
            else
                ψc2 .= ψ.(coils, r2, z2; Bp_fac)
            end

            A[k, :] .= ψc1 .- ψc2

            # weighting
            w = sqrt(cp.weight)
            b[k] *= w
            A[k, :] .*= w

            # store values
            r1_old, z1_old, r2_old, z2_old = r1, z1, r2, z2
            ψc1_old .= ψc1
            ψc2_old .= ψc2
        end
    end

end