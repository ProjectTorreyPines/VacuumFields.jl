module TEQUILA_Ext

import VacuumFields, TEQUILA, MillerExtendedHarmonic, MXHEquilibrium

"""
    encircling_coils(shot::TEQUILA.Shot, n_coils::Int)

Return a Vector of `n_coils` `PointCoil`s distributed outside of `shot`'s boundary
"""
function VacuumFields.encircling_coils(shot::TEQUILA.Shot, n_coils::Int)
    bnd_r, bnd_z = MillerExtendedHarmonic.MXH(shot.surfaces[:, end])()
    return VacuumFields.encircling_coils(bnd_r, bnd_z, n_coils)
end

"""
    boundary_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)

Return a Vector of FluxControlPoint, each with target `ψbound`, at `Npts` equally distributed `fraction_inside` percent inside the the boundary of `shot`
"""
function VacuumFields.boundary_control_points(shot::TEQUILA.Shot, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)
    bnd = MillerExtendedHarmonic.MXH(shot, fraction_inside)
    θs = LinRange(0, 2π, Npts + 1)
    ψtarget = ψbound + TEQUILA.psi_ρθ(shot, fraction_inside, 0.0)
    return [VacuumFields.FluxControlPoint(bnd(θ)..., ψtarget, 1.0 / Npts) for θ in θs[1:end-1]]
end

"""
    boundary_iso_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999; Npts::Integer=99)

Return a Vector of IsoControlPoints, at `Npts` equally distributed `fraction_inside` percent inside the the boundary of `shot`
"""
function VacuumFields.boundary_iso_control_points(shot::TEQUILA.Shot, fraction_inside::Float64=0.999; Npts::Integer=99)
    bnd = MillerExtendedHarmonic.MXH(shot, fraction_inside)
    θs = LinRange(0, 2π, Npts + 1)
    r = [bnd(θ)[1] for θ in θs[1:end-1]]
    z = [bnd(θ)[2] for θ in θs[1:end-1]]
    return VacuumFields.IsoControlPoints(r, z)
end

VacuumFields.plasma_boundary_psi_w_fallback(shot::TEQUILA.Shot, args...) = MXHEquilibrium.Boundary(MXHEquilibrium.MXH(shot, 1.0)()...), 0.0

function VacuumFields.Image(shot::TEQUILA.Shot)
    bnd = @views MillerExtendedHarmonic.MXH(shot.surfaces[:, end])
    return VacuumFields.Image(shot, bnd)
end

function VacuumFields.Image(shot::TEQUILA.Shot, bnd::MillerExtendedHarmonic.MXH; Nb::Integer=100 * (length(bnd.c) + 1))
    image = VacuumFields.Image(shot, bnd(Nb; adaptive=false)...)
    return image
end

end