module VacuumFields

import MillerExtendedHarmonic
import MXHEquilibrium
using Base.Threads
using LinearAlgebra
using RecipesBase
import PlotUtils: cgrad
import FastGaussQuadrature: gausslegendre
import StaticArrays: SMatrix, SVector
import IMAS
import FusionMaterials
import Interpolations
import PlotUtils
using LoopVectorization

const μ₀ = 4e-7 * π
const inv4π = 0.25 / π
const deg2rad = π / 180.0
const default_order = 10

const IMAScoil = IMAS.pf_active__coil
const IMASloop = IMAS.pf_passive__loop
const IMASelement = Union{IMAS.pf_active__coil___element, IMAS.pf_passive__loop___element}
const IMASoutline = Union{IMAS.pf_active__coil___element___geometry__outline, IMAS.pf_passive__loop___element___geometry__outline, NamedTuple}

include("coils.jl")
export AbstractCoil, PointCoil, ParallelogramCoil, QuadCoil, DistributedCoil, MultiCoil, encircling_coils

include("circuit.jl")
export AbstractCircuit, SeriesCircuit, update_coil_currents!

include("integration.jl")

include("image.jl")

include("flux.jl")

include("elliptic.jl")

include("green.jl")

include("control.jl")
export AbstractControlPoint, FluxControlPoint, SaddleControlPoint, IsoControlPoint, FieldControlPoint
export FluxControlPoints, SaddleControlPoints, IsoControlPoints, IsoRefControlPoints, FieldControlPoints
export find_coil_currents!, boundary_control_points

include("fixed2free.jl")
export fixed2free, optimal_λ_regularize

include("mutual.jl")
export mutual, dM_dZ, d2M_dZ2, stability_margin, normalized_growth_rate

include("imas_coils.jl")
export GS_IMAS_pf_active__coil, IMAS_pf_active__coils

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
