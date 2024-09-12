# Generalized functions for specific coil types

function _gfunc(Gfunc::Function, C::PointCoil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=default_order, yorder::Int=default_order)
    return Gfunc(C.R, C.Z, R, Z, scale_factor)
end

function _gfunc(Gfunc::Function, C::Union{ParallelogramCoil, QuadCoil}, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=default_order, yorder::Int=default_order)
    return integrate((X, Y) -> Gfunc(X, Y, R, Z, scale_factor), C; xorder, yorder) / area(C)
end

function _gfunc(Gfunc::Function, C::DistributedCoil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=default_order, yorder::Int=default_order)
    return sum(Gfunc(C.R[k], C.Z[k], R, Z, scale_factor) for k in eachindex(C.R)) / length(C.R)
end

function _gfunc(Gfunc::Function, C::IMAScoil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=default_order, yorder::Int=default_order)
    return sum(_gfunc(Gfunc, element, R, Z, scale_factor; xorder, yorder) for element in elements(C))
end

function _gfunc(Gfunc::Function, element::IMASelement, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=default_order, yorder::Int=default_order)
    return _gfunc(Gfunc, IMAS.outline(element), R, Z, scale_factor; xorder, yorder)
end

function _gfunc(Gfunc::Function, ol::IMASoutline, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=default_order, yorder::Int=default_order)
    return integrate((X, Y) -> Gfunc(X, Y, R, Z, scale_factor), ol; xorder, yorder) / area(ol)
end

function _gfunc(Gfunc::Function, mcoil::MultiCoil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    Ic = current(mcoil)
    Ic == 0.0 && set_current!(mcoil, 1.0)
    G = sum(current(coil) * _gfunc(Gfunc, coil, R, Z, scale_factor; kwargs...) for coil in mcoil.coils) / current(mcoil)
    Ic == 0.0 && set_current!(mcoil, Ic)
    return G
end

# Generalized wrapper functions for all coil types
function Green(coil::Union{AbstractCoil, IMAScoil, IMASelement}, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(Green, coil, R, Z, scale_factor; kwargs...)
end

function dG_dR(coil::Union{AbstractCoil, IMAScoil, IMASelement}, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dG_dR, coil, R, Z, scale_factor; kwargs...)
end

function dG_dZ(coil::Union{AbstractCoil, IMAScoil, IMASelement}, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dG_dZ, coil, R, Z, scale_factor; kwargs...)
end

function d2G_dZ2(coil::Union{AbstractCoil, IMAScoil, IMASelement}, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(d2G_dZ2, coil, R, Z, scale_factor; kwargs...)
end

# Function for circuits
@inline function _gfunc(Pfunc, series::SeriesCircuit, R::Real, Z::Real, scale_factor::Real=1.0;
    COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), Bp_fac::Float64=COCOS.sigma_Bp * (2π)^COCOS.exp_Bp)
    Ic = series.current_per_turn
    try
        update_coil_currents!(series, 1.0)
        return scale_factor * Pfunc(series, R, Z; COCOS, Bp_fac) / (μ₀ * Bp_fac)
    finally
        update_coil_currents!(series, Ic)
    end
end

function Green(series::SeriesCircuit, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(ψ, series, R, Z, scale_factor; kwargs...)
end

function dG_dR(series::SeriesCircuit, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dψ_dR, series, R, Z, scale_factor; kwargs...)
end

function dG_dZ(series::SeriesCircuit, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(dψ_dZ, series, R, Z, scale_factor; kwargs...)
end

function d2G_dZ2(series::SeriesCircuit, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(d2ψ_dZ2, series, R, Z, scale_factor; kwargs...)
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
    return scale_factor * _Green(D, m, Km, Em)
end

_Green(D::Real, m::Real, Km::Real, Em::Real) = inv4π * (2.0 * Em - (2.0 - m) * Km) * sqrt(D)


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
    dG += inv4π * mR * _α(m, Km, dKm, dEm) * sqrt(D) * scale_factor
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
    G = _Green(D, m, Km, Em)
    return scale_factor * _dG_dZ(G, D, DZ_D, m, mZ, Km, dKm, dEm)
end

function _dG_dZ(G::Real, D::Real, DZ_D::Real, m::Real, mZ::Real, Km::Real, dKm::Real, dEm::Real)
    GZ = 0.5 * DZ_D * G
    GZ += inv4π * mZ * _α(m, Km, dKm, dEm) * sqrt(D)
    return GZ
end

# elliptic integral functions used in dG_dR, dG_dZ & d2G_dZ2
function α(X::Real, Y::Real, R::Real, Z::Real)
    _, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    dKm = D_ellipk(m, Km, Em)
    dEm = D_ellipe(m, Km, Em)
    return _α(m, Km, dKm, dEm)
end
_α(m, Km, dKm, dEm) = 2 * dEm - (2.0 - m) * dKm + Km


# Second derivative of Green(X, Y, R, Z) with respect to Z
function d2G_dZ2(X::Real, Y::Real, R::Real, Z::Real, scale_factor::Real=1.0)
    D, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    dKm = D_ellipk(m, Km, Em)
    dEm = D_ellipe(m, Km, Em)
    d2Km = D2_ellipk(m, Km, Em, dKm, dEm)
    d2Em = D2_ellipe(m, Km, Em, dKm, dEm)

    invD = 1.0 / D
    DZ = 2 * (Z - Y)
    DZ_D = DZ * invD
    mZ = -m * DZ_D

    DZZ = 2.0
    mZZ = (m * invD) * (2.0 * DZ * DZ_D - DZZ)
    α = _α(m, Km, dKm, dEm)
    dα_dZ = _dα_dZ(m, mZ, dKm, d2Km, d2Em)

    G = _Green(D, m, Km, Em)
    GZ = _dG_dZ(G, D, DZ_D, m, mZ, Km, dKm, dEm)

    GZZ = 0.5 * invD *  ((DZZ - (DZ * DZ_D)) * G + DZ * GZ)
    GZZ += inv4π * sqrt(D) * ((mZZ + 0.5 * mZ * DZ_D) * α + mZ * dα_dZ)

    return scale_factor * GZZ
end

# elliptic integral functions used in d2G_dZ2
function dα_dZ(X::Real, Y::Real, R::Real, Z::Real)
    D, m = D_m(X, Y, R, Z)
    Km, Em = ellipke(m)
    dKm = D_ellipk(m, Km, Em)
    dEm = D_ellipe(m, Km, Em)
    d2Km = D2_ellipk(m, Km, Em, dKm, dEm)
    d2Em = D2_ellipe(m, Km, Em, dKm, dEm)
    mZ = -2 * m * (Z - Y) / D
    return _dα_dZ(m, mZ, dKm, d2Km, d2Em)
end
_dα_dZ(m, mZ, dKm, d2Km, d2Em) =  mZ * (2.0 * (d2Em + dKm) - (2.0 - m) * d2Km)