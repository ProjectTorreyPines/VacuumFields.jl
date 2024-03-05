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

const coils_D3D_matrix = [
    [8.6080e-01 1.6830e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6140e-01 5.0810e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6280e-01 8.4910e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6110e-01 1.1899e+00 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [1.0041e+00 1.5169e+00 1.3920e-01 1.1940e-01 4.5000e+01 9.0000e+01]
    [2.6124e+00 4.3760e-01 1.7320e-01 1.9460e-01 0.0000e+00 9.2400e+01]
    [2.3733e+00 1.1171e+00 1.8800e-01 1.6920e-01 0.0000e+00 1.0806e+02]
    [1.2518e+00 1.6019e+00 2.3490e-01 8.5100e-02 0.0000e+00 9.0000e+01]
    [1.6890e+00 1.5874e+00 1.6940e-01 1.3310e-01 0.0000e+00 9.0000e+01]
    [8.6080e-01 -1.7370e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6070e-01 -5.1350e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6110e-01 -8.5430e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6300e-01 -1.1957e+00 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [1.0025e+00 -1.5169e+00 1.3920e-01 1.1940e-01 -4.5000e+01 9.0000e+01]
    [2.6124e+00 -4.3760e-01 1.7320e-01 1.9460e-01 0.0000e+00 -9.2400e+01]
    [2.3834e+00 -1.1171e+00 1.8800e-01 1.6920e-01 0.0000e+00 -1.0806e+02]
    [1.2524e+00 -1.6027e+00 2.3490e-01 8.5100e-02 0.0000e+00 9.0000e+01]
    [1.6889e+00 -1.5780e+00 1.6940e-01 1.3310e-01 0.0000e+00 9.0000e+01]]

const coils_D3D_points = [PointCoil(coils_D3D_matrix[i, 1], coils_D3D_matrix[i, 2]) for i in 1:size(coils_D3D_matrix)[1]]
const coils_D3D = [ParallelogramCoil(coils_D3D_matrix[i, :]...) for i in 1:size(coils_D3D_matrix)[1]]

export coils_D3D_points, coils_D3D

end
