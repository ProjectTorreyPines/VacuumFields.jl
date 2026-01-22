using VacuumFields
using Test
using MXHEquilibrium
import EFIT: readg
import DelimitedFiles
using Plots
using MAT

const Ic = 1e5
const resistance = 1e-6
const turns = 10
const Icpt = Ic / turns

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

const coils = [QuadCoil(ParallelogramCoil(pg..., Icpt; resistance, turns)) for pg in eachrow(coils_D3D_matrix)]


function load_par(params)
    R = params[2]
    Z = params[1]
    ΔR = params[4]
    ΔZ = params[3]
    θ₁ = params[5]
    #θ₂ = params[6] - 90.0
    # ??? What's the equation here?
    θ₂ = (params[6] == 0.0) ? 90.0 : params[6]
    return R, Z, ΔR, ΔZ, θ₁, θ₂
end

# Include test files
include("runtests_basic.jl")
include("runtests_ellipke.jl")
include("runtests_inplace.jl")
include("runtests_coil_cache.jl")
