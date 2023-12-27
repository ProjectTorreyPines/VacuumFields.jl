#==============================================================
Define inner products
===============================================================#

function gl_preallocate(N_gl::Integer)
    gξ = zeros(N_gl, N_gl)
    gw = zeros(N_gl, N_gl)
    for i in 1:N_gl
        gξ[1:i, i], gw[1:i, i] = gausslegendre(i)
    end
    return SMatrix{N_gl, N_gl}(gξ),  SMatrix{N_gl, N_gl}(gw)
end

const N_gl = 50
const gξ_pa, gw_pa = gl_preallocate(N_gl)

function integrate(f, lims, order::Integer)
    @assert order <= N_gl
    I = 0.0
    dxdξ = 0.5*(lims[2] - lims[1])
    xavg = 0.5*(lims[2] + lims[1])
    for k in 1:order
        I += f(dxdξ * gξ_pa[k, order] + xavg) * gw_pa[k, order] * dxdξ
    end
    return I
end

function integrate(f, xlims, ylims; xorder=3, yorder=3)
    @assert xorder <= N_gl
    @assert yorder <= N_gl

    dxdξ = 0.5 * (xlims[2] - xlims[1])
    xavg = 0.5 * (xlims[2] + xlims[1])

    dydξ = 0.5 * (ylims[2] - ylims[1])
    yavg = 0.5 * (ylims[2] + ylims[1])

    I = zero(xavg + yavg)

    for i in 1:xorder
        r = dxdξ * gξ_pa[i, xorder] + xavg
        wr = gw_pa[i, xorder] * dxdξ
        Iz = zero(I)
        for j in 1:yorder
            z = dydξ * gξ_pa[j, yorder] + yavg
            wz = gw_pa[j, yorder] * dydξ
            Iz += f(r, z) * wz
        end
        I += Iz * wr
    end
    return I
end

# parameterized R coordinate of parallelogram
function Rplgm(x, y, C::ParallelogramCoil)
    α2 = tan(deg2rad * (C.θ2 + 90.0))
    return Rplgm(x, y, C.R, C.ΔR, C.ΔZ, α2)
end
@inline function Rplgm(x, y, R, ΔR, ΔZ, α2)
    return R + 0.5 * (x * ΔR - y * α2 * ΔZ)
end

# parameterized Z coordinate of parallelogram
function Zplgm(x, y, C::ParallelogramCoil)
    α1 = tan(deg2rad * C.θ1)
    return Zplgm(x, y, C.Z, C.ΔR, C.ΔZ, α1)
end
@inline function Zplgm(x, y, Z, ΔR, ΔZ, α1)
    return Z + 0.5 * (y * ΔZ + x * α1 * ΔR)
end

# integrate over a parallelogram
function integrate(f, C::ParallelogramCoil; xorder=3, yorder=3)
    return integrate(f, C.R, C.Z, C.ΔR, C.ΔZ, C.θ1, C.θ2; xorder, yorder)
end

function integrate(f, R, Z, ΔR, ΔZ, θ1, θ2; xorder=3, yorder=3)
    @assert xorder <= N_gl
    @assert yorder <= N_gl

    α1 = tan(deg2rad * θ1)
    α2 = tan(deg2rad * (θ2 + 90.0))

    @views xs  = gξ_pa[:, xorder]
    @views wxs = gw_pa[:, xorder]
    @views ys  = gξ_pa[:, yorder]
    @views wys = gw_pa[:, yorder]

    g = (x, y) -> f(Rplgm(x, y, R, ΔR, ΔZ, α2), Zplgm(x, y, Z, ΔR, ΔZ, α1))

    J = 0.25 * (1.0 + α1 * α2) * ΔR * ΔZ
    I = J * sum(g(xs[i], ys[j]) * wxs[i] * wys[j] for i in 1:xorder, j in 1:yorder)
    return I
end