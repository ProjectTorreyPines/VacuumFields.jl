using AD_GS
using Test
using EFIT
using Equilibrium

function get_coils()
    coils_D3D = [[ 8.6080e-01  1.6830e-01  5.0800e-02  3.2110e-01  0.0000e+00  9.0000e+01]
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

    Rc = coils_D3D[:,1]
    Zc = coils_D3D[:,2]
    return Rc, Zc
end

const fixed_geqdsk = (@__DIR__)*"/g163303.03170_CHEASE"
const free_geqdsk = (@__DIR__)*"/g163303.03170_EFIT"

@testset "coil_current" begin
    gfixed = transform_cocos(readg(fixed_geqdsk), 5, 1)
    Gfixed = efit(gfixed, clockwise_phi=false)
    
    Rc, Zc = get_coils()
    coils = zip(Rc,Zc)
    currents = fixed_eq_currents(Gfixed,coils)
    println(currents*1e-3)

    gfree = transform_cocos(readg(free_geqdsk), 5, 1)
    Gfree = efit(gfree,clockwise_phi=false)

    check_fixed_eq_currents(Gfixed,coils,currents,Gfree)
end