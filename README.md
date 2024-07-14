# VacuumFields.jl

Closed-boundary solvers do not provide a poloidal field solution in the vacuum region, outside of the last closed flux surface. VacuumFields.jl uses the Green’s function method to match the plasma current contribution with the contributions of external poloidal field coils at the last closed flux surface to determine the homogeneous solution. The total solution is then extended into the vacuum region to get a realistic vacuum solution.

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/VacuumFields.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/VacuumFields.jl/actions/workflows/make_docs.yml/badge.svg)
