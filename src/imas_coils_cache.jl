#= ==================================== =#
#  Cache Management for IMAS Coils      #
#= ==================================== =#
# This file contains cache validation and update logic for GS_IMAS_pf_active__coil.
# Included from imas_coils.jl after struct definitions.

# ===== Element Cache Management =====

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

# ===== Current Cache Management =====

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
