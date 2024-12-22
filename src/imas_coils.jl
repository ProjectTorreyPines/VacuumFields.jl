#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#
mutable struct GS_IMAS_pf_active__coil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    imas::IMAS.pf_active__coil{T1}
    tech::IMAS.build__pf_active__technology{T1}
    time0::Float64
    green_model::Symbol
end

function GS_IMAS_pf_active__coil(
    pfcoil::IMAS.pf_active__coil{T},
    oh_pf_coil_tech::Union{IMAS.build__oh__technology{T},IMAS.build__pf_active__technology{T}},
    green_model::Symbol,
    default_resistance::Float64=1e-6) where {T<:Real}

    # convert everything to build__pf_active__technology so that the `coil_tech`
    # type in GS_IMAS_pf_active__coil is defined at compile time
    coil_tech = IMAS.build__pf_active__technology{T}()
    for field in keys(oh_pf_coil_tech)
        if !isempty(oh_pf_coil_tech, field)
            setproperty!(coil_tech, field, getproperty(oh_pf_coil_tech, field))
        end
    end

    coil = GS_IMAS_pf_active__coil{T,T,T,T}(
        pfcoil,
        coil_tech,
        IMAS.global_time(pfcoil),
        green_model)

    if !ismissing(coil_tech, :material)
        mat_pf = FusionMaterials.Material(coil_tech)
        sigma = mat_pf.electrical_conductivity
        if ismissing(mat_pf) || ismissing(sigma)
            coil.resistance = default_resistance
        else
            coil.resistance = resistance(coil.imas, 1.0 / sigma(; temperature=0.0), :parallel)
        end
    end

    return coil
end

function IMAS_pf_active__coils(dd::IMAS.dd{D}; green_model::Symbol=:quad, zero_currents::Bool=false) where {D<:Real}
    coils = GS_IMAS_pf_active__coil{D,D}[]
    for coil in dd.pf_active.coil
        if zero_currents
            IMAS.@ddtime(coil.current.data = 0.0)   # zero currents for all coils
        end
        if :shaping ∉ (IMAS.index_2_name(coil.function)[f.index] for f in coil.function)
            continue
        end
        if IMAS.is_ohmic_coil(coil)
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        imas_pf_active__coil = GS_IMAS_pf_active__coil(coil, coil_tech, green_model)
        push!(coils, imas_pf_active__coil)
    end
    return coils
end

function imas(coil::GS_IMAS_pf_active__coil)
    return getfield(coil, :imas)
end

function Base.getproperty(coil::GS_IMAS_pf_active__coil{T}, field::Symbol) where {T<:Real}
    pfcoil = getfield(coil, :imas)
    if field ∈ (:r, :z, :width, :height)
        value = getfield(pfcoil.element[1].geometry.rectangle, field)
    elseif field == :current
        # IMAS uses current per turn, GS_IMAS_pf_active__coil uses total current
        value = IMAS.get_time_array(pfcoil.current, :data, getfield(coil, :time0)) .* getproperty(pfcoil.element[1], :turns_with_sign, 1.0)
    elseif field == :resistance
        value = getfield(pfcoil, field)
    elseif field == :turns
        value = abs(getproperty(pfcoil.element[1], :turns_with_sign, 1.0))
    else
        value = getfield(coil, field)
    end
    return value
end

function Base.setproperty!(coil::GS_IMAS_pf_active__coil, field::Symbol, value::Real)
    pfcoil = getfield(coil, :imas)
    if field == :current
        # IMAS uses current per turn, GS_IMAS_pf_active__coil uses total current
        value = value ./ getproperty(pfcoil.element[1], :turns_with_sign, 1.0)
        return IMAS.set_time_array(pfcoil.current, :data, getfield(coil, :time0), value)
    elseif field ∈ (:r, :z, :width, :height)
        return setfield!(pfcoil.element[1].geometry.rectangle, field, value)
    elseif field == :resistance
        setfield!(pfcoil, field, value)
    elseif field == :turns
        val = abs(value) * sign(getproperty(pfcoil.element[1], :turns_with_sign, 1.0))
        setfield!(pfcoil.element[1], :turns_with_sign, val)
    else
        setfield!(coil, field, value)
    end
    return value
end

"""
    Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(Green, coil, R, Z)
end

function dG_dR(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dG_dR, coil, R, Z)
end

function dG_dZ(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dG_dZ, coil, R, Z)
end

function _gfunc(Gfunc::F1, coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=3, yorder::Int=3) where {F1<:Function}
    green_model = getfield(coil, :green_model)

    if green_model == :point # low-fidelity
        oute = IMAS.outline(coil.imas.element[1])
        rc0, zc0 = IMAS.centroid(oute.r, oute.z)
        return Gfunc(rc0, zc0, R, Z)

    elseif green_model == :quad # high-fidelity
        return Gfunc(coil.imas, R, Z)

    else
        error("$(typeof(coil)) green_model is `$(green_model)` but it can only be `:point` or `:quad`")

    end
end


# Mutual inductance

function mutual(C1::GS_IMAS_pf_active__coil, C2::GS_IMAS_pf_active__coil; xorder::Int=3, yorder::Int=3)

    gm1 = getfield(C1, :green_model)
    gm2 = getfield(C2, :green_model)
    @assert gm1 in (:point, :quad)
    @assert gm2 in (:point, :quad)

    if gm1 === :quad && gm2 === :quad
        mutual(C1.imas, C2.imas; xorder, yorder)
    else
        fac = -2π * μ₀ * C1.turns * C2.turns
        if gm1 === :point && gm2 === :point
            return fac * Green(C1.r, C1.z, C2.r, C2.z)
        elseif gm1 === :point
            return fac * Green(C2.imas, C1.r, C1.r)
        else
            return fac * Green(C1.imas, C2.r, C2.r)
        end
    end
end

function mutual(C1::AbstractCoil, C2::GS_IMAS_pf_active__coil; xorder::Int=3, yorder::Int=3)

    green_model = getfield(C2, :green_model)
    if green_model == :point # fastest
        fac = -2π * μ₀ * turns(C1) * C2.turns
        return fac * Green(C1, C2.r, C2.z)

    elseif green_model == :quad # high-fidelity
        return mutual(C1, C2.imas; xorder, yorder)

    else
        error("$(typeof(C2)) green_model can only be (in order of accuracy) :quad and :point")
    end
end


function _pfunc(Pfunc, image::Image, C::GS_IMAS_pf_active__coil, δZ;
    COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), xorder::Int=3, yorder::Int=3)

    green_model = getfield(C, :green_model)
    if green_model == :point # fastest
        PC = PointCoil(C.r, C.z; C.turns)
        return _pfunc(Pfunc, image, PC, δZ; COCOS)

    elseif green_model == :quad # high-fidelity
        return _pfunc(Pfunc, image, C.imas, δZ; COCOS, xorder, yorder)

    else
        error("$(typeof(C)) green_model can only be (in order of accuracy) :quad and :point")
    end
end

@recipe function plot_coil(coil::GS_IMAS_pf_active__coil)
    @series begin
        #seriestype := :scatter
        #marker --> :circle
        label --> ""
        #[coil.r], [coil.z]
        coil.imas
    end
end

@recipe function plot_coil(coils::AbstractVector{<:GS_IMAS_pf_active__coil})
    for (k, coil) in enumerate(coils)
        @series begin
            primary := (k == 1)
            aspect_ratio := :equal
            coil
        end
    end
end