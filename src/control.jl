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
        seriestype := :scatter
        aspect_ratio := :equal
        label --> ""
        markersize := cp.weight * 10
        [cp.R], [cp.Z]
    end
end

@recipe function plot_FluxControlPoint(cp::FluxControlPoint)
    @series begin
        marker := :circle
        seriestype := :scatter
        aspect_ratio := :equal
        label --> ""
        markersize := cp.weight * 10
        [cp.R], [cp.Z]
    end
end

@recipe function plot_IsoControlPoint(cp::IsoControlPoint)
    @series begin
        marker := :diamond
        linewidth := 0.1
        aspect_ratio := :equal
        label --> ""
        markersize := cp.weight * 10
        [cp.R1, cp.R2], [cp.Z1, cp.Z2]
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
function IsoControlPoints(Rs::AbstractVector{T}, Zs::AbstractVector{T}) where {T<:Real}
    @assert length(Zs) == length(Rs)
    is_closed = sqrt((Rs[1] - Rs[end])^2 + (Zs[1] - Zs[end])^2) < 1E-6
    if is_closed
        Rs = @view Rs[1:end-1]
        Zs = @view Zs[1:end-1]
    end
    N = length(Rs)
    i_max = argmax(Rs)
    i_min = argmin(Rs)
    iso_cps = IsoControlPoint{T}[]
    weight = 1.0 / (N - 1)
    # connect i_min to i_max
    push!(iso_cps, IsoControlPoint{T}(Rs[i_min], Zs[i_min], Rs[i_max], Zs[i_max], weight))
    # then interlace the rest
    for k in (i_max+1):(i_max+N-1) # exclude i_max explicitly
        ck = IMAS.index_circular(N, k)
        (ck == i_min) && continue # skip i_min... we did it first
        if mod(k, 2) == mod(i_max, 2)
            push!(iso_cps, IsoControlPoint{T}(Rs[ck], Zs[ck], Rs[i_min], Zs[i_min], weight))
        else
            push!(iso_cps, IsoControlPoint{T}(Rs[ck], Zs[ck], Rs[i_max], Zs[i_max], weight))
        end
    end
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
    boundary_iso_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999; Npts::Integer=99)

Return a Vector of IsoControlPoints, at `Npts` equally distributed `fraction_inside` percent inside the the boundary of `EQfixed`
"""
function boundary_iso_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999; Npts::Integer=99)
    ψ0, ψb = MXHEquilibrium.psi_limits(EQfixed)
    Sp = MXHEquilibrium.flux_surface(EQfixed, fraction_inside * (ψb - ψ0) + ψ0; n_interp=Npts)
    return VacuumFields.IsoControlPoints(Sp.r, Sp.z)
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
        Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1],
        cocos=MXHEquilibrium.cocos(EQ),
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

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
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1],
    cocos=MXHEquilibrium.cocos(EQ),
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    return find_coil_currents!(coils, EQ, Image(EQ); flux_cps, saddle_cps, iso_cps, ψbound, fixed_coils, λ_regularize, Sb, A)
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
        Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1],
        cocos=MXHEquilibrium.cocos(EQ),
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

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
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1],
    cocos=MXHEquilibrium.cocos(EQ),
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    b = init_b(EQ, image; flux_cps, saddle_cps, iso_cps, ψbound, Sb)
    offset_b!(b; flux_cps, saddle_cps, iso_cps, fixed_coils, cocos)
    return _find_coil_currents!(coils, A, b; λ_regularize)
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        ψpl::Interpolations.AbstractInterpolation,
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        ψbound::Real=0.0,
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        cocos=MXHEquilibrium.cocos(11),
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Assumes flux from plasma current given by the ψpl interpolation and a `ψbound` flux at the boundary `Sb`

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    ψpl::Interpolations.AbstractInterpolation,
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11),
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    dψpl_dR = (r, z) -> Interpolations.gradient(ψpl, r, z)[1]
    dψpl_dZ = (r, z) -> Interpolations.gradient(ψpl, r, z)[2]
    return find_coil_currents!(coils, ψpl, dψpl_dR, dψpl_dZ;
                               flux_cps, saddle_cps, iso_cps,
                               fixed_coils, λ_regularize, cocos, A)

end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil},
        ψpl::Union{Function,Interpolations.AbstractInterpolation},
        dψpl_dR::Union{Function,Interpolations.AbstractInterpolation},
        dψpl_dZ::Union{Function,Interpolations.AbstractInterpolation};
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        cocos=MXHEquilibrium.cocos(11),
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

Find the currents for `coils` that best match (leas-squares) the control points provided by `flux_cps`, `saddle_cps`, and `iso_cps`

Assumes flux and R/Z gradients from plasma current given by the ψpl, dψpl_dR, and dψpl_dZ interpolations and a `ψbound` flux at the boundary `Sb`

Optionally assumes flux from additional `fixed_coils`, whose currents will not change

`λ_regularize` provides regularization in the least-squares fitting
"""
function find_coil_currents!(
    coils::AbstractVector{<:AbstractCoil},
    ψpl::Union{Function,Interpolations.AbstractInterpolation},
    dψpl_dR::Union{Function,Interpolations.AbstractInterpolation},
    dψpl_dZ::Union{Function,Interpolations.AbstractInterpolation};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    λ_regularize::Float64=0.0,
    cocos=MXHEquilibrium.cocos(11),
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    b = init_b(ψpl, dψpl_dR, dψpl_dZ; flux_cps, saddle_cps, iso_cps)
    offset_b!(b; flux_cps, saddle_cps, iso_cps, fixed_coils, cocos)
    _find_coil_currents!(coils, A, b; λ_regularize)
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
        Sb=nothing,
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

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
    Sb=nothing,
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    return find_coil_currents!(coils; flux_cps, saddle_cps, iso_cps, fixed_coils, λ_regularize, cocos, A)
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
        Sb=nothing,
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

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
    Sb=nothing,
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    return find_coil_currents!(coils; flux_cps, saddle_cps, iso_cps, fixed_coils, λ_regularize, cocos, A)
end

"""
    find_coil_currents!(
        coils::AbstractVector{<:AbstractCoil};
        flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
        saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
        iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
        fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
        λ_regularize::Float64=0.0,
        cocos=MXHEquilibrium.cocos(11),
        A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

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
    cocos=MXHEquilibrium.cocos(11),
    A::AbstractMatrix{<:Real}=define_A(coils; flux_cps, saddle_cps, iso_cps, cocos))

    N = length(flux_cps) + 2 * length(saddle_cps) + length(iso_cps)
    b = zeros(N)
    offset_b!(b; flux_cps, saddle_cps, iso_cps, fixed_coils, cocos)
    return _find_coil_currents!(coils, A, b; λ_regularize)
end

function _find_coil_currents!(coils::AbstractVector{<:AbstractCoil},
                              A::AbstractMatrix{<:Real},
                              b::AbstractVector{<:Real};
                              λ_regularize::Float64=0.0)

    @assert size(A) == (length(b), length(coils))

    # Least-squares solve for coil currents
    if λ_regularize > 0
        # Least-squares with regularization
        # https://www.youtube.com/watch?v=9BckeGN0sF0
        Ic0 = reg_solve(A, b, λ_regularize / length(coils)^2)
    else
        Ic0 = A \ b
    end
    #Ic0 here is "total current" in each coil
    cost = norm(A * Ic0 .- b) / norm(b)

    # update values of coils current
    for (k, coil) in enumerate(coils)
        Ic0[k] /= turns(coil) # convert to current per turn
        set_current_per_turn!(coil, Ic0[k])
    end

    return Ic0, cost
end

function init_b(
    EQ::MXHEquilibrium.AbstractEquilibrium,
    image::Image;
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    ψbound::Real=0.0,
    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])

    _, ψb = MXHEquilibrium.psi_limits(EQ)

    ψpl = (r, z) -> (MXHEquilibrium.in_boundary(Sb, (r, z)) ? EQ(r, z) : ψb) + ψbound - ψb - ψ(image, r, z)
    dψpl_dR = (r, z) -> -dψ_dR(image, r, z)
    dψpl_dZ = (r, z) -> -dψ_dZ(image, r, z)

    return init_b(ψpl, dψpl_dR, dψpl_dZ; flux_cps, saddle_cps, iso_cps)

end

function init_b(
    ψpl::Union{Function,Interpolations.AbstractInterpolation},
    dψpl_dR::Union{Function,Interpolations.AbstractInterpolation},
    dψpl_dZ::Union{Function,Interpolations.AbstractInterpolation};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[])

    Nflux = length(flux_cps)
    Nsaddle = length(saddle_cps)
    Niso = length(iso_cps)
    N = Nflux + 2 * Nsaddle + Niso

    T = eltype(iso_cps).parameters[1]
    b = zeros(T, N)

    for i in eachindex(flux_cps)
        cp = flux_cps[i]
        r = cp.R
        z = cp.Z

        # subtract plasma current contribution
        b[i] -= ψpl(r, z)
    end

    for i in eachindex(saddle_cps)
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
    for i in eachindex(iso_cps)
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

function offset_b!(b::AbstractVector{T};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],
    cocos=MXHEquilibrium.cocos(11)) where {T<:Real}

    Nflux = length(flux_cps)
    Nsaddle = length(saddle_cps)

    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp

    for i in eachindex(flux_cps)
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

        # weighting
        w = sqrt(cp.weight)
        b[i] *= w
    end

    for i in eachindex(saddle_cps)
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

        # weighting
        w = sqrt(cp.weight)
        b[ir:iz] .*= w
    end

    if !isempty(iso_cps)
        S = promote_type(T, eltype(iso_cps).parameters[1])
        r1_old, z1_old, r2_old, z2_old = -one(S), zero(S), -one(S), zero(S)
        ψf1_old, ψf2_old = zero(S), zero(S)
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

            # weighting
            w = sqrt(cp.weight)
            b[k] *= w

        end
    end

end

function define_A(coils::AbstractVector{<:AbstractCoil};
    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],
    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],
    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],
    cocos=MXHEquilibrium.cocos(11))

    Nflux = length(flux_cps)
    Nsaddle = length(saddle_cps)
    Niso = length(iso_cps)
    N = Nflux + 2 * Nsaddle + Niso

    A = zeros(N, length(coils))

    Bp_fac = cocos.sigma_Bp * (2π)^cocos.exp_Bp

    # First reset current in coils to unity
    for coil in coils
        set_current_per_turn!(coil, 1.0/turns(coil))
    end

    for i in eachindex(flux_cps)
        cp = flux_cps[i]
        r = cp.R
        z = cp.Z

        # RHS

        # Build matrix relating coil Green's functions to boundary points
        A[i, :] .= ψ.(coils, r, z; Bp_fac)

        # weighting
        w = sqrt(cp.weight)
        A[i, :] .*= w
    end

    for i in eachindex(saddle_cps)
        cp = saddle_cps[i]
        r = cp.R
        z = cp.Z

        ir = Nflux + 2i - 1
        iz = Nflux + 2i

        # Build matrix relating coil Green's functions to boundary points
        A[ir, :] .= dψ_dR.(coils, r, z; Bp_fac)
        A[iz, :] .= dψ_dZ.(coils, r, z; Bp_fac)

        # weighting
        w = sqrt(cp.weight)
        A[ir:iz, :] .*= w
    end

    if !isempty(iso_cps)
        T = eltype(iso_cps).parameters[1]
        r1_old, z1_old, r2_old, z2_old = -one(T), zero(T), -one(T), zero(T)
        Ncoil = length(coils)
        ψc1_old, ψc2_old = Vector{T}(undef, Ncoil), Vector{T}(undef, Ncoil)
        ψc1, ψc2 = Vector{T}(undef, Ncoil), Vector{T}(undef, Ncoil)
        for i in eachindex(iso_cps)
            cp = iso_cps[i]
            r1, z1 = cp.R1, cp.Z1
            r2, z2 = cp.R2, cp.Z2

            k = Nflux + 2Nsaddle + i

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
            A[k, :] .*= w

            # store values
            r1_old, z1_old, r2_old, z2_old = r1, z1, r2, z2
            ψc1_old .= ψc1
            ψc2_old .= ψc2
        end
    end

    return A
end