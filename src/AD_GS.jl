__precompile__()

using Equilibrium
using Plots
using Trapz
using Unzip
using SpecialFunctions
using Optim

module AD_GS

const μ₀ = 4e-7*π

export fixed_eq_currents, check_fixed_eq_currents

include("coil_currents.jl")

end
