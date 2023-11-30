const inv4π = 0.25 / π

@inline function D_m(X::Real, Y::Real, R::Real, Z::Real)
    D = (X + R)^2 + (Y - Z)^2
    m = 4.0 * X * R / D # this is k^2
    return D, m
end

function Green(C::ParallelogramCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return Green(DistributedCoil(C), R, Z, scale_factor)
end

function Green(C::DistributedCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return sum(Green(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end

function Green(C::PointCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return Green(C.R, C.Z, R, Z, scale_factor)
end

@inline function Green(X::Real, Y::Real, R::Real, Z::Real, scale_factor::Real=1.0)
    D, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    return inv4π * (2.0 * Em - (2.0 - m) * Km) * sqrt(D) * scale_factor
end


# Derivative of Green(X, Y, R, Z) with respect to R

function dG_dR(C::ParallelogramCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return dG_dR(DistributedCoil(C), R, Z, scale_factor)
end

function dG_dR(C::DistributedCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return sum(dG_dR(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end

function dG_dR(C::PointCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return dG_dR(C.R, C.Z, R, Z, scale_factor)
end

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

function dG_dZ(C::ParallelogramCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return dG_dZ(DistributedCoil(C), R, Z, scale_factor)
end

function dG_dZ(C::DistributedCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return sum(dG_dZ(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end

function dG_dZ(C::PointCoil, R::Real, Z::Real, scale_factor::Real=1.0)
    return dG_dZ(C.R, C.Z, R, Z, scale_factor)
end

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