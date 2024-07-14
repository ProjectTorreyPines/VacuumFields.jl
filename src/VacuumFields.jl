module VacuumFields

import MillerExtendedHarmonic
import MXHEquilibrium
using Base.Threads
using LinearAlgebra
using RecipesBase
import PlotUtils: cgrad
import TEQUILA
import FastGaussQuadrature: gausslegendre
import StaticArrays: SMatrix, SVector
import IMAS

const μ₀ = 4e-7 * π
const inv4π = 0.25 / π
const deg2rad = π / 180.0
const default_order = 10

const IMAScoil = IMAS.pf_active__coil
const IMASelement = IMAS.pf_active__coil___element
const IMASoutline = Union{IMAS.pf_active__coil___element___geometry__outline, NamedTuple}

include("coils.jl")

include("integration.jl")

include("image.jl")

include("flux.jl")

include("elliptic.jl")

include("green.jl")

include("control.jl")

include("fixed2free.jl")

include("mutual.jl")

export AbstractCoil, PointCoil, DistributedCoil, coil, encircling_coils, plot_coil
export AbstractControlPoint, FluxControlPoint, SaddleControlPoint, find_coil_currents!
export fixed2free, optimal_λ_regularize, encircling_fixed2free
export mutual, dM_dZ, d2M_dZ2, stability_margin, normalized_growth_rate

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
