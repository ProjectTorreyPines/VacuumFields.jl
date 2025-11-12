#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#

"""
    ElementCache{T}

Cached derived quantities from IMAS coil data to avoid repeated allocation and computation.
Each element in a coil has one ElementCache.
Includes source geometry for cache validation.
Mutable to enable in-place updates for zero-allocation cache invalidation.
"""
mutable struct ElementCache{T<:Real}
    turns_with_sign::T
    outline_r::Vector{T}
    outline_z::Vector{T}
    source_geometry::NamedTuple{(:r, :z, :width, :height), NTuple{4, T}}
end

"""
    CurrentCache{T}

Cached current_per_turn value for a specific time point.
Stores the time when current was evaluated and its corresponding current_per_turn.
"""
struct CurrentCache{T<:Real}
    time::T
    current_per_turn::T
end

mutable struct GS_IMAS_pf_active__coil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    imas::IMAS.pf_active__coil{T1}
    tech::IMAS.build__pf_active__technology{T1}
    time0::Float64
    green_model::Symbol
    _elements_cache::Vector{ElementCache{T1}}
    _current_cache::CurrentCache{T1}

    # Inner constructor (default initialization)
    function GS_IMAS_pf_active__coil{T1,T2,T3,T4}(
        imas::IMAS.pf_active__coil{T1},
        tech::IMAS.build__pf_active__technology{T1},
        time0::Float64,
        green_model::Symbol) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}
        empty_cache = ElementCache{T1}[]
        invalid_current = CurrentCache{T1}(T1(NaN), T1(NaN))
        new{T1,T2,T3,T4}(imas, tech, time0, green_model, empty_cache, invalid_current)
    end

    # Inner constructor (with element cache)
    function GS_IMAS_pf_active__coil{T1,T2,T3,T4}(
        imas::IMAS.pf_active__coil{T1},
        tech::IMAS.build__pf_active__technology{T1},
        time0::Float64,
        green_model::Symbol,
        elements_cache::Vector{ElementCache{T1}}) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}
        invalid_current = CurrentCache{T1}(T1(NaN), T1(NaN))
        new{T1,T2,T3,T4}(imas, tech, time0, green_model, elements_cache, invalid_current)
    end

    # Inner constructor (with element cache and current cache)
    function GS_IMAS_pf_active__coil{T1,T2,T3,T4}(
        imas::IMAS.pf_active__coil{T1},
        tech::IMAS.build__pf_active__technology{T1},
        time0::Float64,
        green_model::Symbol,
        elements_cache::Vector{ElementCache{T1}},
        current_cache::CurrentCache{T1}) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}
        new{T1,T2,T3,T4}(imas, tech, time0, green_model, elements_cache, current_cache)
    end
end

# Helper functions for element cache management

"""
    is_cache_valid(elem_cache::ElementCache, element) -> Bool

Check if a single ElementCache is still valid for the given element.
Returns false if geometry or turns_with_sign changed.
"""
@inline function is_cache_valid(elem_cache::ElementCache, element)
    rect = element.geometry.rectangle

    # Compare geometry fields
    current_geom = (r=rect.r, z=rect.z, width=rect.width, height=rect.height)
    current_geom != elem_cache.source_geometry && return false

    # turns_with_sign comparison
    element.turns_with_sign != elem_cache.turns_with_sign && return false

    return true
end

"""
    is_cache_valid(elements_cache::Vector{ElementCache{T}}, elements) -> Bool

Check if all cached element data is still valid by comparing source geometry.
Returns false if cache is empty, element count changed, or any geometry/turns changed.
"""
@inline function is_cache_valid(elements_cache::Vector{ElementCache{T}}, elements) where {T}
    isempty(elements_cache) && return false
    length(elements_cache) != length(elements) && return false

    for (k, element) in enumerate(elements)
        !is_cache_valid(elements_cache[k], element) && return false
    end

    return true
end

"""Create ElementCache from a single element (internal helper)"""
@inline function _create_element_cache(element, ::Type{T}) where {T}
    rect = element.geometry.rectangle
    ol = IMAS.outline(element)
    return ElementCache{T}(
        element.turns_with_sign,
        ol.r,
        ol.z,
        (r=rect.r, z=rect.z, width=rect.width, height=rect.height)
    )
end

"""
    update_cache!(cache::ElementCache, element)

Update ElementCache in-place from element data. Minimizes allocation by reusing vectors.
"""
@inline function update_cache!(cache::ElementCache{T}, element) where {T}
    rect = element.geometry.rectangle
    ol = IMAS.outline(element)

    # Update scalar
    cache.turns_with_sign = element.turns_with_sign

    # Update vectors in-place (resize + copy)
    resize!(cache.outline_r, length(ol.r))
    copyto!(cache.outline_r, ol.r)

    resize!(cache.outline_z, length(ol.z))
    copyto!(cache.outline_z, ol.z)

    # Update geometry (small NamedTuple, cheap to reassign)
    cache.source_geometry = (r=rect.r, z=rect.z, width=rect.width, height=rect.height)

    return cache
end

"""
    ensure_valid_elements_cache!(coil)

Ensure all elements have valid cache. Updates only invalid elements in-place.
If cache is empty or size mismatched, regenerates entire cache.
This is the main entry point for cache management.
"""
function ensure_valid_elements_cache!(coil::GS_IMAS_pf_active__coil{T}) where {T}
    # If cache is empty or size mismatched, regenerate entire cache
    if isempty(coil._elements_cache) || length(coil._elements_cache) != length(coil.imas.element)
        # Full regeneration
        setfield!(coil, :_elements_cache, [_create_element_cache(el, T) for el in coil.imas.element])
        return coil
    end

    # Otherwise, validate and update only invalid elements IN-PLACE (zero allocation)
    for (k, element) in enumerate(coil.imas.element)
        if !is_cache_valid(coil._elements_cache[k], element)
            update_cache!(coil._elements_cache[k], element)
        end
    end

    return coil
end

"""Clear cached element data"""
function clear_elements_cache!(coil::GS_IMAS_pf_active__coil)
    setfield!(coil, :_elements_cache, ElementCache{eltype(coil._elements_cache)}[])
    return coil
end

# Helper functions for current cache management

"""
    is_cache_valid(current_cache::CurrentCache, time) -> Bool

Check if cached current is valid for the specified time.
Returns true if cache is valid (not NaN) and time matches, false otherwise.
"""
@inline function is_cache_valid(current_cache::CurrentCache, time::Real)
    return !isnan(current_cache.time) && current_cache.time == time
end


"""
    cache_current!(coil, time, current_per_turn)

Cache current_per_turn value for the specified time.
"""
function cache_current!(coil::GS_IMAS_pf_active__coil{T}, time::Real, current_per_turn::Real) where {T}
    setfield!(coil, :_current_cache, CurrentCache{T}(T(time), T(current_per_turn)))
    return coil
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
    elseif field == :current_per_turn
        time0 = getfield(coil, :time0)
        if is_cache_valid(coil._current_cache, time0)
            # Cache hit - use cached value
            value = coil._current_cache.current_per_turn
        else
            # Cache miss - fetch from IMAS and cache it
            value = IMAS.get_time_array(pfcoil.current, :data, time0)
            cache_current!(coil, time0, value)
        end
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
    if field == :current_per_turn
        time0 = getfield(coil, :time0)
        # Update cache to reflect the new value
        cache_current!(coil, time0, value)
        return IMAS.set_time_array(pfcoil.current, :data, time0, value)
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
    return _dispatch_green(Green, coil, R, Z)
end

function dG_dR(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _dispatch_green(dG_dR, coil, R, Z)
end

function dG_dZ(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _dispatch_green(dG_dZ, coil, R, Z)
end

function _dispatch_green(Gfunc::F1, coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=3, yorder::Int=3) where {F1<:Function}
    green_model = getfield(coil, :green_model)

    if green_model == :point # low-fidelity
        oute = IMAS.outline(coil.imas.element[1])
        rc0, zc0 = IMAS.centroid(oute.r, oute.z)
        return Gfunc(rc0, zc0, R, Z, scale_factor) * coil.turns

    elseif green_model == :quad # high-fidelity
        # Ensure cache is valid (updates only invalid elements)
        ensure_valid_elements_cache!(coil)

        # Compute using cached data
        sum_val = 0.0
        for (k, element) in enumerate(coil.imas.element)
            sum_val += _gfunc(Gfunc, element, R, Z,
                            coil._elements_cache[k].outline_r, coil._elements_cache[k].outline_z, coil._elements_cache[k].turns_with_sign,
                            scale_factor; xorder, yorder)
        end

        return sum_val
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
        fac = -2π * μ₀
        if gm1 === :point && gm2 === :point
            return fac * C1.turns * C2.turns * Green(C1.r, C1.z, C2.r, C2.z)
        elseif gm1 === :point
            # C2.turns inside Green
            return fac * C1.turns * Green(C2.imas, C1.r, C1.r; xorder, yorder)
        else
            # C1.turns inside Green
            return fac * C2.turns * Green(C1.imas, C2.r, C2.r; xorder, yorder)
        end
    end
end

function mutual(C1::AbstractCoil, C2::GS_IMAS_pf_active__coil; xorder::Int=3, yorder::Int=3)

    green_model = getfield(C2, :green_model)
    if green_model == :point # fastest
        # C1.turns inside Green
        fac = -2π * μ₀ * C2.turns
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