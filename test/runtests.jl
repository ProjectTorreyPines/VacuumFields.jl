using AD_GS
using Test
using EFIT
using Equilibrium

const fixed_geqdsk = (@__DIR__)*"/g163303.03170_CHEASE"
const free_geqdsk = (@__DIR__)*"/g163303.03170_EFIT"
const good_currents = [-197968.56551716285, -175313.42410039817, -52907.88075466734, -256837.87396872052, 485346.0233782668, -214060.41838723494, -151469.54405769982, -330872.8548105091, -488.24739294769034, -194334.88010685574, -136597.64936162502, -159445.93012418432, 16601.967623453867, 53881.80298066415, -237456.52831624495, -230790.53129979773, 85387.28458689603, -23680.644329152925]

@testset "coil_current" begin
    gfixed = transform_cocos(readg(fixed_geqdsk), 5, 1)
    Gfixed = efit(gfixed, clockwise_phi=false)
    
    currents = fixed_eq_currents(Gfixed,coils_D3D)
    @test currents ≈ good_currents

    gfree = transform_cocos(readg(free_geqdsk), 5, 1)
    Gfree = efit(gfree,clockwise_phi=false)

    check_fixed_eq_currents(Gfixed,coils,currents,Gfree)
end
