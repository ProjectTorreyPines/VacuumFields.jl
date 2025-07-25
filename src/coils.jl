"""
    AbstractCoil

Abstract coil type
"""
abstract type AbstractCoil end

"""
    AbstractSingleCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real}

Abstract coil type
"""
abstract type AbstractSingleCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractCoil end
(::Type{T})(R, Z; resistance=0.0, turns=1) where {T<:AbstractSingleCoil} = T(R, Z, 0.0, resistance, turns)
(::Type{T})(R, Z, current_per_turn; resistance=0.0, turns=1) where {T<:AbstractSingleCoil} = T(R, Z, current_per_turn, resistance, turns)

"""
    PointCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}

Point filament coil at scalar (R, Z) location
"""
mutable struct PointCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::T1
    Z::T1
    current_per_turn::T2
    resistance::T3
    turns::T4
end

"""
    ParallelogramCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}

Parallelogram coil with the R, Z, ΔR, ΔZ, θ1, θ2 formalism (as used by EFIT, for example).
Here θ1 and θ2 are the shear angles along the x- and y-axes, respectively, in degrees.
"""
mutable struct ParallelogramCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::T1
    Z::T1
    ΔR::T1
    ΔZ::T1
    θ1::T1
    θ2::T1
    current_per_turn::T2
    resistance::T3
    turns::T4
end

"""
    QuadCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}

Quadrilateral coil with counter-clockwise corners (starting from lower left) at R and Z
"""
mutable struct QuadCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real,VT1<:AbstractVector{T1}} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::VT1
    Z::VT1
    current_per_turn::T2
    resistance::T3
    turns::T4

    function QuadCoil(R::VT1, Z::VT1, current_per_turn::T2, resistance::T3, turns::T4) where {VT1,T2,T3,T4}
        @assert length(R) == length(Z) == 4
        points = reverse!(IMAS.convex_hull!(collect(zip(R, Z)); closed_polygon=false)) # reverse to make ccw
        a, b = zip(points...)
        R = VT1(collect(a))
        Z = VT1(collect(b))
        return new{eltype(VT1),T2,T3,T4,VT1}(R, Z, current_per_turn, resistance, turns)
    end
end

"""
    DistributedCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}

Coil consisting of distributed filaments
"""
mutable struct DistributedCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::Vector{T1}
    Z::Vector{T1}
    current_per_turn::T2
    resistance::T3
    turns::T4
end

"""
    MultiCoil{SC<:AbstractSingleCoil} <: AbstractCoil

A coil consisting of multiple coils linke in series
"""
mutable struct MultiCoil{SC<:AbstractSingleCoil} <: AbstractCoil
    coils::Vector{SC}
    orientation::Vector{Int}
    function MultiCoil(coils::Vector{SC}, orientation::Vector{Int}; copy_coils::Bool=true) where {SC<:AbstractSingleCoil}
        @assert length(coils) == length(orientation)
        total_turns = sum(turns(coil) for coil in coils)
        Ic = sum(current_per_turn(coil) * turns(coil) * orientation[k] for (k, coil) in enumerate(coils))
        new_coils = copy_coils ? deepcopy(coils) : coils
        mcoil = new{SC}(new_coils, orientation)
        set_current_per_turn!(mcoil, Ic / total_turns)
        return mcoil
    end
end

current_per_turn(coil::AbstractSingleCoil) = coil.current_per_turn

elements(coil::Union{IMAScoil,IMASloop}) = Iterators.filter(!isempty, coil.element)

# BCL 2/27/24
# N.B.: Not sure about sign with turns and such
function current_per_turn(coil::IMAScoil)
    if ismissing(coil.current, :data)
        return 0.0
    end
    return IMAS.@ddtime(coil.current.data)
end

function current_per_turn(loop::IMASloop)
    if ismissing(loop, :current)
        return 0.0
    end
    return IMAS.@ddtime(loop.current)
end

@inline function set_current_per_turn!(coil::AbstractSingleCoil, current_per_turn::Real)
    return coil.current_per_turn = current_per_turn
end

function set_current_per_turn!(coil::IMAScoil, current_per_turn::Real)
    IMAS.@ddtime(coil.current.data = current_per_turn)
    return coil
end

resistance(coil::AbstractSingleCoil) = coil.resistance

function resistance(coil::IMAScoil{T}) where {T<:Real}
    return ismissing(coil, :resistance) ? zero(T) : coil.resistance
end

turns(coil::AbstractSingleCoil) = coil.turns
# VacuumFields turns are sign-less
function turns(coil::IMAScoil; with_sign=false)
    Nt = sum(turns(element) for element in elements(coil))
    return with_sign ? Nt : abs(Nt)
end

function turns(element::IMASelement; with_sign=false)
    T = eltype(element)
    isempty(element) && return zero(T)
    Nt = ismissing(element, :turns_with_sign) ? one(T) : element.turns_with_sign
    return with_sign ? Nt : abs(Nt)
end

"""
    PointCoil(R, Z, current_per_turn=0.0; resistance=0.0, turns=1)

Return PointCoil, a point filament coil at scalar (R, Z) location
"""
PointCoil(R, Z,  current_per_turn=0.0; resistance=0.0, turns=1) = PointCoil(R, Z, current_per_turn, resistance, turns)



"""
    ParallelogramCoil(R, Z, ΔR, ΔZ, θ1, θ2, current_per_turn=0.0; resistance=0.0, turns=1)

Construct a ParallelogramCoil
"""
function ParallelogramCoil(R, Z, ΔR, ΔZ, θ1, θ2, current_per_turn=0.0; resistance=0.0, turns=1)
    R, Z, ΔR, ΔZ, θ1, θ2 = promote(R, Z, ΔR, ΔZ, θ1, θ2)
    return ParallelogramCoil(R, Z, ΔR, ΔZ, θ1, θ2, current_per_turn, resistance, turns)
end

area(C::ParallelogramCoil) = area(C.ΔR, C.ΔZ, C.θ1, C.θ2)

function area(ΔR, ΔZ, θ1, θ2)
    α1 = tan(deg2rad * θ1)
    α2 = tan(deg2rad * (θ2 + 90.0))
    return (1.0 + α1 * α2) * ΔR * ΔZ
end


"""
    QuadCoil(R, Z, current_per_turn=0.0; resistance=0.0, turns=1)

Construct a QuadCoil from corners defined by `R` and `Z`
"""
QuadCoil(R, Z, current_per_turn=0.0; resistance=0.0, turns=1) = QuadCoil(R, Z, current_per_turn, resistance, turns)

"""
    QuadCoil(pc::ParallelogramCoil)

Construct a QuadCoil from an existing ParallelogramCoil
"""
function QuadCoil(pc::ParallelogramCoil)
    x = SVector(-1.0, 1.0, 1.0, -1.0)
    y = SVector(-1.0, -1.0, 1.0, 1.0)
    R = VacuumFields.Rplgm.(x, y, Ref(pc))
    Z = VacuumFields.Zplgm.(x, y, Ref(pc))
    return QuadCoil(R, Z, current_per_turn(pc), resistance(pc), turns(pc))
end

function QuadCoil(elm::IMASelement, current_per_turn::Real=0.0, resistance_per_turn::Real=0.0)
    R, Z = IMAS.outline(elm)
    @assert length(R) == length(Z) == 4
    Nt = turns(elm)
    return QuadCoil(R, Z, current_per_turn, resistance_per_turn * Nt, Nt)
end

function area(R::AbstractVector{<:Real}, Z::AbstractVector{<:Real})
    @assert length(R) == length(Z) == 4
    A = R[1] * Z[2] + R[2] * Z[3] + R[3] * Z[4] + R[4] * Z[1]
    A -= R[2] * Z[1] + R[3] * Z[2] + R[4] * Z[3] + R[1] * Z[4]
    return 0.5 * abs(A)
end

area(C::QuadCoil) = area(C.R, C.Z)
area(coil::IMAScoil) = sum(area(element) for element in elements(coil))
area(element::IMASelement) = area(IMAS.outline(element))
area(ol::IMASoutline) = area(ol.r, ol.z)

# compute the resistance given a resistitivity
function resistance(coil::Union{ParallelogramCoil,QuadCoil,IMASelement}, resistivity::Real)
    return 2π * turns(coil)^2 * resistivity / integrate((R, Z) -> 1.0 / R, coil)
end

function resistance(coil::IMAScoil, resistivity::Real, element_connection::Symbol=:series)
    if element_connection === :series
        eta = sum(resistance(element, resistivity) for element in elements(coil))
    elseif element_connection === :parallel
        eta = 1.0 / sum(1.0 / resistance(element, resistivity) for element in elements(coil))
    else
        error("element_connection should be :series or :parallel")
    end
    return eta
end

function resistance(mcoil::MultiCoil, resistivity::Real)
    return sum(resistance(coil, resistivity) for coil in mcoil.coils)
end

function set_resistance!(coil::Union{ParallelogramCoil,QuadCoil}, resistivity::Real)
    return coil.resistance = resistance(coil, resistivity)
end

function set_resistance!(coil::IMAScoil, resistivity::Real, element_connection::Symbol=:series)
    return coil.resistance = resistance(coil, resistivity, element_connection)
end

function set_resistance!(mcoil::MultiCoil, resistivity::Real)
    for coil in mcoil.coils
        coil.resistance = resistance(coil, resistivity)
    end
end


# compute the resistivity of a coil based on it's resistance
function resistivity(coil::Union{ParallelogramCoil,QuadCoil,IMASelement})
    return resistance(coil) * integrate((R, Z) -> 1.0 / R, coil) / (2π * turns(coil)^2)
end

function resistivity(coil::IMAScoil, element_connection::Symbol)
    if element_connection === :series
        rho = resistance(coil, element_connection) / sum(resistance(element, 1.0) for element in elements(coil))
    elseif element_connection === :parallel
        rho = resistance(coil, element_connection) * sum(1.0 / resistance(element, 1.0) for element in elements(coil))
    else
        error("element_connection should be :series or :parallel")
    end
    return rho
end

function resistivity(mcoil::MultiCoil)
    return resistance(mcoil) / sum(resistance(coil, 1.0) for coil in mcoil.coils)
end


"""
    DistributedCoil(R::Vector{<:Real}, Z::Vector{<:Real}, current_per_turn=0.0; resistance=0.0, turns=1)

Construct a DistributedCoil from filaments defined by `R` and `Z` vectors
"""
DistributedCoil(R::Vector{<:Real}, Z::Vector{<:Real}, current_per_turn=0.0; resistance=0.0, turns=1) = DistributedCoil(R, Z, current_per_turn, resistance, turns)

"""
    DistributedParallelogramCoil(Rc::T1, Zc::T1, ΔR::T1, ΔZ::T1, θ1::T1, θ2::T1, current_per_turn::Real=0.0; spacing::Real=0.01, turns::Real=1) where {T1<:Real}

Create a DistributedCoil of filaments with `spacing` separation within
a parallelogram defined by the R, Z, ΔR, ΔZ, θ1, θ2 formalism (as used by EFIT, for example)

NOTE: if spacing <= 0.0 then current filaments are placed at the vertices
"""
function DistributedParallelogramCoil(Rc::T1, Zc::T1, ΔR::T1, ΔZ::T1, θ1::T1, θ2::T1, current_per_turn::Real=0.0; spacing::Real=0.01, turns::Real=1) where {T1<:Real}
    if spacing <= 0.0
        dR = [-0.5 * ΔR, 0.5 * ΔR]
        dZ = [-0.5 * ΔZ, 0.5 * ΔZ]
    else
        dR = LinRange(-0.5 * ΔR, 0.5 * ΔR, round(Int, 1.0 + ΔR / spacing, RoundUp))
        dZ = LinRange(-0.5 * ΔZ, 0.5 * ΔZ, round(Int, 1.0 + ΔZ / spacing, RoundUp))
    end

    α1 = tan(π * θ1 / 180.0)
    α2 = tan(π * (θ2 + 90.0) / 180.0)

    Nr = length(dR)
    Nz = length(dZ)
    N = Nr * Nz
    Z = Array{Float64,1}(undef, N)
    R = Array{Float64,1}(undef, N)
    for i in 1:Nr
        for j in 1:Nz
            k = j + (i - 1) * Nz
            @inbounds R[k] = Rc + dR[i] - α2 * dZ[j]
            @inbounds Z[k] = Zc + dZ[j] + α1 * dR[i]
        end
    end
    return DistributedCoil(R, Z, current_per_turn; turns)
end

"""
    DistributedCoil(C::ParallelogramCoil; spacing::Real=0.01)

Create a DistributedCoil from an existing ParallelogramCoil, with filaments of `spacing` separation

NOTE: if spacing <= 0.0 then current filaments are placed at the vertices of parallelogram
"""
function DistributedCoil(C::ParallelogramCoil; spacing::Real=0.01)
    return DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ1, C.θ2, current_per_turn(C); spacing, turns=turns(C))
end


# MultiCoil

function MultiCoil(icoil::IMAScoil)
    total_turns = turns(icoil)
    resistance_per_turn = resistance(icoil) / total_turns

    Icpt = current_per_turn(icoil)
    orientation = [Int(sign(turns(elm; with_sign=true))) for elm in elements(icoil)]
    coils = [QuadCoil(elm, Icpt * orientation[k], resistance_per_turn) for (k, elm) in enumerate(elements(icoil))]
    return MultiCoil(coils, orientation)
end

function MultiCoil(loop::IMASloop)
    total_turns = sum(abs(ismissing(element, :turns_with_sign) ? 1.0 : element.turns_with_sign) for element in elements(loop))
    if !ismissing(loop, :resistance)
        resistance_per_turn = loop.resistance / total_turns
    else
        resistance_per_turn = 0.0
    end

    # I'm assuming that pf_passive is like pf_active and loop.current is current/turn like coil.current is
    Icpt = current_per_turn(loop)
    orientation = [Int(sign(turns(elm; with_sign=true))) for elm in elements(loop)]
    coils = [QuadCoil(elm, Icpt * orientation[k], resistance_per_turn) for (k, elm) in enumerate(elements(loop))]
    if resistance_per_turn == 0.0 && !ismissing(loop, :resistivity)
        for coil in coils
            coil.resistance = resistance(coil, loop.resistivity)
        end
    end

    return MultiCoil(coils, orientation)
end

function MultiCoils(dd::IMAS.dd{D}; load_pf_active::Bool=true, active_only::Bool=false, load_pf_passive::Bool=false) where {D<:Real}
    if load_pf_active
        active_coils = MultiCoils(dd.pf_active; active_only)
        !load_pf_passive && return active_coils
    end
    if load_pf_passive
        passive_coils = MultiCoils(dd.pf_passive)
        !load_pf_active && return passive_coils
    end
    isempty(active_coils) && return passive_coils
    isempty(passive_coils) && return active_coils
    return vcat(active_coils, passive_coils)
end

function MultiCoils(pf_active::IMAS.pf_active; active_only::Bool=false)
    pf_coils = active_only ? filter(is_active, pf_active.coil) : pf_active.coil
    return [MultiCoil(coil) for coil in pf_coils]
end

function MultiCoils(pf_passive::IMAS.pf_passive)
    return [MultiCoil(loop) for loop in pf_passive.loop]
end

function current_per_turn(mcoil::MultiCoil)
    Icpt = current_per_turn(mcoil.coils[1]) * mcoil.orientation[1]
    for (k, coil) in enumerate(mcoil.coils)
        @assert current_per_turn(coil) * mcoil.orientation[k] == Icpt "All coils in MultiCoil must have the same current per turn"
    end
    return Icpt
end
resistance(mcoil::MultiCoil) = sum(resistance(coil) for coil in mcoil.coils)
turns(mcoil::MultiCoil) = sum(turns(coil) for coil in mcoil.coils)

function set_current_per_turn!(mcoil::MultiCoil, current_per_turn::Real)
    for (k, coil) in enumerate(mcoil.coils)
        set_current_per_turn!(coil, current_per_turn * mcoil.orientation[k])
    end
    return mcoil
end

function update_orientation!(mcoil::MultiCoil, orientation::Vector{Int})
    @assert length(mcoil.coils) == length(orientation)
    Icpt = current_per_turn(mcoil)
    mcoil.orientation = orientation
    set_current_per_turn!(mcoil, Icpt)
    return mcoil
end

# IMAS related functions

function is_active(coil::IMAScoil)
    funcs = (IMAS.index_2_name(coil.function)[f.index] for f in coil.function)
    return :shaping in funcs
end

function update_currents!(icoils::AbstractVector{<:IMAScoil}, coils::AbstractVector{<:AbstractCoil}; active_only::Bool=false)
    if active_only
        @assert sum(is_active, icoils) == length(coils)
        k = 0
        for icoil in icoils
            if is_active(icoil)
                k += 1
                VacuumFields.set_current_per_turn!(icoil, VacuumFields.current_per_turn(coils[k]))
            end
        end
    else
        @assert length(icoils) == length(coils)
        for (k, icoil) in enumerate(icoils)
            VacuumFields.set_current_per_turn!(icoil, VacuumFields.current_per_turn(coils[k]))
        end
    end
    return
end

"""
    encircling_coils(EQfixed::MXHEquilibrium.AbstractEquilibrium, n_coils::Int)

Return a Vector of `n_coils` `PointCoil`s distributed outside of `EQfixed`'s boundary
"""
function encircling_coils(EQfixed::MXHEquilibrium.AbstractEquilibrium, n_coils::Int)
    bnd = MXHEquilibrium.plasma_boundary(EQfixed; precision=0.0)
    if bnd === nothing
        # if the original boundary specified in EQfixed does not close, then find LCFS boundary
        bnd = MXHEquilibrium.plasma_boundary(EQfixed)
    end

    return encircling_coils(bnd.r, bnd.z, n_coils)
end

"""
    encircling_coils(bnd_r::AbstractVector{T}, bnd_z::AbstractVector{T}, n_coils::Int) where {T<:Real}

Return a Vector of `n_coils` `PointCoil`s distributed outside of closed boundary defined by `bnd_r` and `bnd_z`
"""
function encircling_coils(bnd_r::AbstractVector{T}, bnd_z::AbstractVector{T}, n_coils::Int) where {T<:Real}
    mxh = MillerExtendedHarmonic.MXH(bnd_r, bnd_z, 2)
    mxh.R0 = mxh.R0 .* (1.0 + mxh.ϵ)
    mxh.ϵ = 0.999
    mxh.κ *= 1.2
    Θ = LinRange(0, 2π, n_coils + 1)[1:end-1]
    return [PointCoil(r, z) for (r, z) in mxh.(Θ)]
end

# ============== #
#   Convex hull  #
# ============== #
function IMAS.convex_hull(C::ParallelogramCoil)
    C = DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ1, C.θ2; spacing=0.0)
    return collect(zip(C.R, C.Z))
end

function IMAS.convex_hull(C::DistributedCoil)
    return IMAS.convex_hull(C.R, C.Z; closed_polygon=true)
end

# ======== #
#   Plot   #
# ======== #
@recipe function plot_coil(C::ParallelogramCoil)
    return DistributedCoil(C)
end

@recipe function plot_coil(C::DistributedCoil)
    hull = IMAS.convex_hull(C)
    R = [r for (r, z) in hull]
    Z = [z for (r, z) in hull]
    @series begin
        seriestype --> :shape
        color --> :black
        alpha --> 0.2
        label --> ""
        R, Z
    end
end

@recipe function plot_coil(C::QuadCoil)
    @series begin
        seriestype --> :shape
        color --> :black
        alpha --> 0.2
        label --> ""
        C.R, C.Z
    end
end

@recipe function plot_coil(C::PointCoil)
    @series begin
        seriestype --> :scatter
        marker --> :circle
        markercolor --> :black
        [C.R], [C.Z]
    end
end

max_current(coil::AbstractCoil) = abs(current_per_turn(coil) * turns(coil))
max_current(mcoil::MultiCoil) = isempty(mcoil.coils) ? 0.0 : maximum(abs, current_per_turn(coil) * turns(coil) for coil in mcoil.coils)
max_current(icoil::IMAScoil) = ismissing(icoil.current, :data) ? 0.0 : (abs(IMAS.@ddtime(icoil.current.data) * maximum(turns(elm) for elm in icoil.element)))

@recipe function plot_coils(coils::AbstractVector{<:AbstractCoil}; color_by=:current, cname=:diverging)
    if !isempty(coils)
        @assert color_by in (nothing, :current)
        if color_by === :current
            alpha --> nothing
            currents = [current_per_turn(coil) * turns(coil) for coil in coils]
            CURRENT = maximum(max_current(coil) for coil in coils)
            if CURRENT > 1e6
                currents = currents .* 1e-6
                CURRENT = CURRENT .* 1e-6
                c_unit = "MA"
            elseif CURRENT > 1e3
                currents = currents .* 1e-3
                CURRENT = CURRENT .* 1e-3
                c_unit = "kA"
            else
                c_unit = "A"
            end

            @series begin
                colorbar_title := " \nPF currents [$c_unit]"
                seriestype --> :scatter
                color --> cname
                clim --> (-CURRENT, CURRENT)
                marker_z --> [-CURRENT, CURRENT]
                [(NaN, NaN), (NaN, NaN)]
            end
        end

        for (k, coil) in enumerate(coils)
            @series begin
                if color_by === :current
                    colorbar_title := " \nPF currents [$c_unit]"
                    colorbar_entry := false
                    if CURRENT == 0.0
                        color --> :black
                    else
                        current_color_index = (currents[k] + CURRENT) / (2 * CURRENT)
                        color --> PlotUtils.cgrad(cname)[current_color_index]
                    end
                else
                    primary := (k == 1)
                end
                aspect_ratio := :equal
                label --> ""
                coil
            end
        end
    end
end

@recipe function plot_coils(mcoil::MultiCoil)
    for (k, coil) in enumerate(mcoil.coils)
        @series begin
            primary := (k == 1)
            aspect_ratio := :equal
            coil
        end
    end
end