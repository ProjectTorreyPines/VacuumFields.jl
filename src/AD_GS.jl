__precompile__()

module AD_GS

using Base: Real, Integer
using Equilibrium
using Plots
using Trapz
using Unzip
import SpecialFunctions
using Optim
using Base.Threads
using LinearAlgebra
import LazySets: convex_hull, VPolygon, Singleton
import Memoize

const μ₀ = 4e-7*π
const inv2π = 1.0/(2π)

include("coil_currents.jl")

include("elliptic.jl")

export plot_coil_flux, check_fixed_eq_currents, fixed2free
export fixed_eq_currents, ψp_on_fixed_eq_boundary, currents_to_match_ψp
export AbstractCoil, PointCoil, DistributedCoil, coil
export plot_coils, plot_coil

const coils_D3D_matrix = [[ 8.6080e-01  1.6830e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 8.6140e-01  5.0810e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 8.6280e-01  8.4910e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 8.6110e-01  1.1899e+00  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 1.0041e+00  1.5169e+00  1.3920e-01  1.1940e-01  4.5000e+01  9.0000e+01]
                          [ 2.6124e+00  4.3760e-01  1.7320e-01  1.9460e-01  0.0000e+00  9.2400e+01]
                          [ 2.3733e+00  1.1171e+00  1.8800e-01  1.6920e-01  0.0000e+00  1.0806e+02]
                          [ 1.2518e+00  1.6019e+00  2.3490e-01  8.5100e-02  0.0000e+00  9.0000e+01]
                          [ 1.6890e+00  1.5874e+00  1.6940e-01  1.3310e-01  0.0000e+00  9.0000e+01]
                          [ 8.6080e-01 -1.7370e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 8.6070e-01 -5.1350e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 8.6110e-01 -8.5430e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 8.6300e-01 -1.1957e+00  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
                          [ 1.0025e+00 -1.5169e+00  1.3920e-01  1.1940e-01 -4.5000e+01  9.0000e+01]
                          [ 2.6124e+00 -4.3760e-01  1.7320e-01  1.9460e-01  0.0000e+00 -9.2400e+01]
                          [ 2.3834e+00 -1.1171e+00  1.8800e-01  1.6920e-01  0.0000e+00 -1.0806e+02]
                          [ 1.2524e+00 -1.6027e+00  2.3490e-01  8.5100e-02  0.0000e+00  9.0000e+01]
                          [ 1.6889e+00 -1.5780e+00  1.6940e-01  1.3310e-01  0.0000e+00  9.0000e+01]]

Ncoil = size(coils_D3D_matrix)[1]

const coils_D3D_points = [PointCoil(coils_D3D_matrix[i,1],coils_D3D_matrix[i,2]) for i in 1:Ncoil]
const coils_D3D = [ParallelogramCoil(coils_D3D_matrix[i,:]...) for i in 1:Ncoil]

export coils_D3D_points, coils_D3D

end
