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
(::Type{T})(R, Z, current; resistance=0.0, turns=1) where {T<:AbstractSingleCoil} = T(R, Z, current, resistance, turns)

"""
    PointCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}

Point filament coil at scalar (R, Z) location
"""
mutable struct PointCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::T1
    Z::T1
    current::T2
    resistance::T3
    turns::T4
end

current(coil::AbstractSingleCoil) = coil.current

elements(coil::Union{IMAScoil, IMASloop}) = Iterators.filter(!isempty, coil.element)

# BCL 2/27/24
# N.B.: Not sure about sign with turns and such
function current(coil::IMAScoil)
    if ismissing(coil.current, :data)
        return 0.0
    end
    cur_per_turn = IMAS.@ddtime(coil.current.data)
    if cur_per_turn == 0.0
        return cur_per_turn
    else
        return cur_per_turn * sum(turns(element; with_sign=true) for element in elements(coil))
    end
end

function current(loop::IMASloop)
    if ismissing(loop, :current)
        return 0.0
    end
    cur_per_turn = IMAS.@ddtime(loop.current)
    if cur_per_turn == 0.0
        return cur_per_turn
    else
        return cur_per_turn * sum(turns(element; with_sign=true) for element in elements(loop))
    end
end

@inline function set_current!(coil::AbstractSingleCoil, current::Real)
    coil.current = current
end

function set_current!(coil::IMAScoil, current::Real)
    turns = sum(element.turns_with_sign for element in elements(coil))
    ipt = current / turns
    IMAS.@ddtime(coil.current.data = ipt)
    return coil
end

resistance(coil::AbstractSingleCoil) = coil.resistance

function resistance(coil::IMAScoil{T}) where {T<:Real}
    try
        return coil.resistance
    catch e
        if e isa IMAS.IMASmissingDataException
            return zero(T)
        else
            rethrow(e)
        end
    end
end

turns(coil::AbstractSingleCoil) = coil.turns
# VacuumFields turns are sign-less
function turns(coil::IMAScoil; with_sign=false)
    Nt = sum(element.turns_with_sign for element in elements(coil))
    return with_sign ? Nt : abs(Nt)
end

function turns(element::IMASelement; with_sign=false)
    isempty(element) && return 0.0
    Nt = ismissing(element, :turns_with_sign) ? 1.0 : element.turns_with_sign
    return with_sign ? Nt : abs(Nt)
end

"""
    PointCoil(R, Z, current=0.0; resistance=0.0, turns=1)

Return PointCoil, a point filament coil at scalar (R, Z) location
"""
PointCoil(R, Z, current=0.0; resistance=0.0, turns=1) = PointCoil(R, Z, current, resistance, turns)

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
    current::T2
    resistance::T3
    turns::T4
end

"""
    ParallelogramCoil(R, Z, ΔR, ΔZ, θ1, θ2, current=0.0; resistance=0.0, turns=1)

Construct a ParallelogramCoil
"""
function ParallelogramCoil(R, Z, ΔR, ΔZ, θ1, θ2, current=0.0; resistance=0.0, turns=1)
    R, Z, ΔR, ΔZ, θ1, θ2 = promote(R, Z, ΔR, ΔZ, θ1, θ2)
    ParallelogramCoil(R, Z, ΔR, ΔZ, θ1, θ2, current, resistance, turns)
end

area(C::ParallelogramCoil) = area(C.ΔR, C.ΔZ, C.θ1, C.θ2)

function area(ΔR, ΔZ, θ1, θ2)
    α1 = tan(deg2rad * θ1)
    α2 = tan(deg2rad * (θ2 + 90.0))
    return (1.0 + α1 * α2) * ΔR * ΔZ
end

"""
    QuadCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}

Quadrilateral coil with counter-clockwise corners (starting from lower left) at R and Z
"""
mutable struct QuadCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real,VT1<:AbstractVector{T1}} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::VT1
    Z::VT1
    current::T2
    resistance::T3
    turns::T4
    function QuadCoil(R::VT1, Z::VT1, current::T2, resistance::T3, turns::T4) where {VT1,T2,T3,T4}
        @assert length(R) == length(Z) == 4
        points = reverse!(grahamscan!(collect(zip(R, Z)))) # reverse to make ccw
        a, b = zip(points...)
        R = VT1(collect(a))
        Z = VT1(collect(b))
        return new{eltype(VT1),T2,T3,T4,VT1}(R, Z, current, resistance, turns)
    end
end

"""
    QuadCoil(R, Z, current=0.0; resistance=0.0, turns=1)

Construct a QuadCoil from corners defined by `R` and `Z`
"""
QuadCoil(R, Z, current=0.0; resistance=0.0, turns=1) = QuadCoil(R, Z, current, resistance, turns)

"""
    QuadCoil(pc::ParallelogramCoil)

Construct a QuadCoil from an existing ParallelogramCoil
"""
function QuadCoil(pc::ParallelogramCoil)
    x = SVector(-1.0, 1.0, 1.0, -1.0)
    y = SVector(-1.0, -1.0, 1.0, 1.0)
    R = VacuumFields.Rplgm.(x, y, Ref(pc))
    Z = VacuumFields.Zplgm.(x, y, Ref(pc))
    return QuadCoil(R, Z, pc.current; pc.resistance, pc.turns)
end


function QuadCoil(elm::IMASelement, current_per_turn::Real = 0.0, resistance_per_turn::Real = 0.0)
    R, Z = IMAS.outline(elm)
    @assert length(R) == length(Z) == 4
    Nt = turns(elm)
    return QuadCoil(R, Z, current_per_turn * Nt, resistance_per_turn * Nt, Nt)
end

function area(R::AbstractVector{<:Real}, Z::AbstractVector{<:Real})
    @assert length(R) == length(Z) == 4
    A = R[1] * Z[2] + R[2] * Z[3] + R[3] * Z[4] + R[4] * Z[1]
    A -= R[2] * Z[1] + R[3] * Z[2] + R[4] * Z[3] + R[1] * Z[4]
    return 0.5 * A
end

area(C::QuadCoil) = area(C.R, C.Z)
area(coil::IMAScoil) = sum(area(element) for element in elements(coil))
area(element::IMASelement) = area(IMAS.outline(element))
area(ol::IMASoutline) = area(ol.r, ol.z)

# compute the resistance given a resistitivity
function resistance(C::Union{ParallelogramCoil,QuadCoil}, resistivity::Real)
    return 2π * C.turns^2 * resistivity / integrate((R, Z) -> 1.0 / R, C)
end

function resistance(coil::IMAScoil, resistivity::Real, element_connection::Symbol)
    if element_connection === :series
        eta = sum(resistance(element, resistivity) for element in elements(coil))
    elseif element_connection === :parallel
        eta = 1.0 / sum(1.0 / resistance(element, resistivity) for element in elements(coil))
    else
        error("element_connection should be :series or :parallel")
    end
    return eta
end

function resistance(element::IMASelement, resistivity::Real)
    return 2π * turns(element)^2 * resistivity / integrate((R, Z) -> 1.0 / R, element)
end

"""
    DistributedCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}

Coil consisting of distributed filaments
"""
mutable struct DistributedCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}
    R::Vector{T1}
    Z::Vector{T1}
    current::T2
    resistance::T3
    turns::T4
end

"""
    DistributedCoil(R::Vector{<:Real}, Z::Vector{<:Real}, current=0.0; resistance=0.0, turns=1)

Construct a DistributedCoil from filaments defined by `R` and `Z` vectors
"""
DistributedCoil(R::Vector{<:Real}, Z::Vector{<:Real}, current=0.0; resistance=0.0, turns=1) = DistributedCoil(R, Z, current, resistance, turns)

"""
    DistributedParallelogramCoil(Rc::T1, Zc::T1, ΔR::T1, ΔZ::T1, θ1::T1, θ2::T1, current::Real=0.0; spacing::Real=0.01, turns::Real=1) where {T1<:Real}

Create a DistributedCoil of filaments with `spacing` separation within
    a parallelogram defined by the R, Z, ΔR, ΔZ, θ1, θ2 formalism (as used by EFIT, for example)

NOTE: if spacing <= 0.0 then current filaments are placed at the vertices
"""
function DistributedParallelogramCoil(Rc::T1, Zc::T1, ΔR::T1, ΔZ::T1, θ1::T1, θ2::T1, current::Real=0.0; spacing::Real=0.01, turns::Real=1) where {T1<:Real}
    if spacing <= 0.0
        dR = [-0.5 * ΔR, 0.5 * ΔR]
        dZ = [-0.5 * ΔZ, 0.5 * ΔZ]
    else
        dR = LinRange(-0.5 * ΔR, 0.5 * ΔR, Int(ceil(1.0 + ΔR / spacing)))
        dZ = LinRange(-0.5 * ΔZ, 0.5 * ΔZ, Int(ceil(1.0 + ΔZ / spacing)))
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
    return DistributedCoil(R, Z, current; turns)
end

"""
    DistributedCoil(C::ParallelogramCoil; spacing::Real=0.01)

Create a DistributedCoil from an existing ParallelogramCoil, with filaments of `spacing` separation

NOTE: if spacing <= 0.0 then current filaments are placed at the vertices of parallelogram
"""
function DistributedCoil(C::ParallelogramCoil; spacing::Real=0.01)
    return DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ1, C.θ2, C.current; spacing, C.turns)
end


mutable struct MultiCoil{SC<:AbstractSingleCoil} <: AbstractCoil
    coils::Vector{SC}
    function MultiCoil(coils::Vector{SC}; copy_coils::Bool=true) where {SC<:AbstractSingleCoil}
        Ic = sum(current(coil) for coil in coils)
        new_coils = copy_coils ? deepcopy(coils) : coils
        mcoil = new{SC}(new_coils)
        set_current!(mcoil, Ic)
        return mcoil
    end
end

function MultiCoil(icoil::IMAScoil)
    total_turns = turns(icoil)
    resistance_per_turn = resistance(icoil) / total_turns

    current_per_turn = current(icoil) / total_turns
    return MultiCoil([QuadCoil(elm, current_per_turn, resistance_per_turn) for elm in elements(icoil)])
end

function MultiCoil(loop::IMASloop)
    total_turns = sum(abs(ismissing(element, :turns_with_sign) ? 1.0 : element.turns_with_sign) for element in elements(loop))
    if !ismissing(loop, :resistance)
        resistance_per_turn = loop.resistance / total_turns
    else
        resistance_per_turn = 0.0
    end

    # I'm assuming that pf_passive is like pf_passive and loop.current is current/turn like coil.current is
    current_per_turn = ismissing(loop, :current) ? 0.0 : IMAS.@ddtime(loop.current)

    coils = [QuadCoil(elm, current_per_turn, resistance_per_turn) for elm in elements(loop)]
    if resistance_per_turn == 0.0 && !ismissing(loop, :resistivity)
        for coil in coils
            coil.resistance = resistance(coil, loop.resistivity)
        end
    end

    return MultiCoil(coils)
end

function is_active(coil::IMAScoil)
    funcs = (IMAS.index_2_name(coil.function)[f.index] for f in coil.function)
    return (:shaping in funcs) && current(coil) != 0.0
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
    return vcat(active_coils, passive_coils)
end

function MultiCoils(pf_active::IMAS.pf_active; active_only::Bool=false)
    pf_coils = active_only ? filter(is_active, pf_active.coil) : pf_active.coil
    return [MultiCoil(coil) for coil in pf_coils]
end

function MultiCoils(pf_passive::IMAS.pf_passive)
    return [MultiCoil(loop) for loop in pf_passive.loop]
end

current(mcoil::MultiCoil) = sum(current(coil) for coil in mcoil.coils)
resistance(mcoil::MultiCoil) = sum(resistance(coil) for coil in mcoil.coils)
turns(mcoil::MultiCoil) = sum(turns(coil) for coil in mcoil.coils)

function set_current!(mcoil::MultiCoil, current::Real)
    total_turns = turns(mcoil)
    current_per_turn = current / total_turns
    for coil in mcoil.coils
        Nt = turns(coil)
        set_current!(coil, current_per_turn * Nt)
    end
    return mcoil
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
    encircling_coils(shot::TEQUILA.Shot, n_coils::Int)

Return a Vector of `n_coils` `PointCoil`s distributed outside of `shot`'s boundary
"""
function encircling_coils(shot::TEQUILA.Shot, n_coils::Int)
    bnd_r, bnd_z = MillerExtendedHarmonic.MXH(shot.surfaces[:, end])()
    return encircling_coils(bnd_r, bnd_z, n_coils)
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
@inline function points_isless(p::AbstractVector{T}, q::AbstractVector{T}) where {T}
    return p[1] < q[1] || (p[1] == q[1] && p[2] < q[2])
end

@inline function points_isless(p::Tuple{T,T}, q::Tuple{T,T}) where {T}
    return p[1] < q[1] || (p[1] == q[1] && p[2] < q[2])
end

@inline function isrightturn(p::AbstractVector{T}, q::AbstractVector{T}, r::AbstractVector{T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) < 0.0
end

@inline function isrightturn(p::Tuple{T,T}, q::Tuple{T,T}, r::Tuple{T,T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) < 0.0
end

function halfhull(points::AbstractVector)
    halfhull = similar(points)
    n = 0
    for p in points
        while n > 1 && !isrightturn(halfhull[n-1], halfhull[n], p)
            n -= 1
        end
        n += 1
        halfhull[n] = p
    end
    return view(halfhull, 1:n)
end

function grahamscan!(points::AbstractVector)
    sort!(points; lt=points_isless)
    upperhull = halfhull(points)
    reverse!(points)
    lowerhull = halfhull(points)
    return [upperhull; lowerhull[2:end-1]]
end

function convex_hull!(xy_points::AbstractVector; closed_polygon::Bool)
    hull = grahamscan!(xy_points)
    if closed_polygon && !isempty(hull)
        return push!(hull, hull[1])
    else
        return hull
    end
end

function convex_hull(xy_points::AbstractVector; closed_polygon::Bool)
    return convex_hull!(deepcopy(xy_points); closed_polygon)
end

function convex_hull(x::AbstractVector{T}, y::AbstractVector{T}; closed_polygon::Bool) where {T}
    xy_points = [(xx, yx) for (xx, yx) in zip(x, y)]
    return convex_hull!(xy_points; closed_polygon)
end

function convex_hull(C::ParallelogramCoil)
    C = DistributedParallelogramCoil(C.R, C.Z, C.ΔR, C.ΔZ, C.θ1, C.θ2; spacing=0.0)
    return collect(zip(C.R, C.Z))
end

function convex_hull(C::DistributedCoil)
    return convex_hull(C.R, C.Z; closed_polygon=true)
end

# ======== #
#   Plot   #
# ======== #
@recipe function plot_coil(C::ParallelogramCoil)
    return DistributedCoil(C)
end

@recipe function plot_coil(C::DistributedCoil)
    hull = convex_hull(C)
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

@recipe function plot_coils(coils::AbstractVector{<:Union{AbstractCoil,IMAScoil}})
    for (k, coil) in enumerate(coils)
        @series begin
            primary := k == 1
            aspect_ratio := :equal
            label --> ""
            coil
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