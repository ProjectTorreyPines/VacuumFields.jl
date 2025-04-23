const CoilVectorType = AbstractVector{<:Union{VacuumFields.AbstractSingleCoil, IMAS.pf_active__coil, IMAS.pf_active__coil___element}}

abstract type AbstractCircuit end

mutable struct SeriesCircuit{CV<:CoilVectorType, T<:Real} <: AbstractCircuit
    coils::CV
    current_per_turn::T
    signs::Vector{Int}
    function SeriesCircuit(coils::CV, current_per_turn::T, signs; copy_coils::Bool=true) where {CV<:CoilVectorType, T<:Real}
        @assert length(coils) == length(signs)
        new_coils = copy_coils ? deepcopy(coils) : coils
        update_coil_currents!(new_coils, current_per_turn, signs)
        return new{CV, T}(new_coils, current_per_turn, signs)
    end
end

function update_coil_currents!(series::SeriesCircuit, current_per_turn::Real=series.current_per_turn)
    series.current_per_turn = current_per_turn
    update_coil_currents!(series.coils, current_per_turn, series.signs)
    return series
end

function update_coil_currents!(coils::CoilVectorType, current_per_turn::Real, signs::Vector{Int})
    for (k, coil) in enumerate(coils)
        set_current_per_turn!(coil, signs[k] * current_per_turn)
    end
    return coils
end

resistance(series::SeriesCircuit) = sum(resistance(coil) for coil in series.coils)