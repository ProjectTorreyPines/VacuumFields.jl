__precompile__()

module AD_GS

using Equilibrium
using Plots
using Trapz
using Unzip
using SpecialFunctions
using Optim

const μ₀ = 4e-7*π

include("coil_currents.jl")
export fixed_eq_currents, check_fixed_eq_currents

end
