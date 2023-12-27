# Generalized functions for specific coil types

function _gfunc(Gfunc::Function, C::PointCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return Gfunc(C.R, C.Z, R, Z, scale_factor)
end

function _gfunc(Gfunc::Function, C::ParallelogramCoil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=3, yorder::Int=3)
    return integrate((X, Y) -> Gfunc(X, Y, R, Z, scale_factor), C; xorder, yorder) / area(C)
end

function _gfunc(Gfunc::Function, C::DistributedCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return sum(Gfunc(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end


# Generalized wrapper functions for all coil types
function Green(coil::AbstractCoil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(Green, coil, R, Z, scale_factor; kwargs...)
end

function dG_dR(coil::AbstractCoil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dG_dR, coil, R, Z, scale_factor; kwargs...)
end

function dG_dZ(coil::AbstractCoil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dG_dZ, coil, R, Z, scale_factor; kwargs...)
end


# Point-to-point Green's functions

@inline function D_m(X::Real, Y::Real, R::Real, Z::Real)
    D = (X + R)^2 + (Y - Z)^2
    m = 4.0 * X * R / D # this is k^2
    return D, m
end

# Green(X, Y, R, Z)
@inline function Green(X::Real, Y::Real, R::Real, Z::Real, scale_factor::Real=1.0)
    D, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    return inv4π * (2.0 * Em - (2.0 - m) * Km) * sqrt(D) * scale_factor
end

# Derivative of Green(X, Y, R, Z) with respect to R
function dG_dR(X::Real, Y::Real, R::Real, Z::Real, scale_factor::Real=1.0)
    D, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    dKm = D_ellipk(m, Km, Em)
    dEm = D_ellipe(m, Km, Em)

    DR = 2 * (X + R)
    DR_D = DR / D
    mR = m * ((1.0 / R) - DR_D)

    dG = 0.5 * DR_D * Green(X, Y, R, Z, scale_factor)
    dG += inv4π * mR * (2 * dEm - (2.0 - m) * dKm + Km) * sqrt(D) * scale_factor
    return dG
end


# Derivative of Green(X, Y, R, Z) with respect to Z
function dG_dZ(X::Real, Y::Real, R::Real, Z::Real, scale_factor::Real=1.0)
    D, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    dKm = D_ellipk(m, Km, Em)
    dEm = D_ellipe(m, Km, Em)

    DZ = 2 * (Z - Y)
    DZ_D = DZ / D
    mZ = -m * DZ_D

    dG = 0.5 * DZ_D * Green(X, Y, R, Z, scale_factor)
    dG += inv4π * mZ * (2 * dEm - (2.0 - m) * dKm + Km) * sqrt(D) * scale_factor
    return dG
end