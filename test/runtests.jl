using AD_GS
using Test
using EFIT
using Equilibrium

const fixed_geqdsk = (@__DIR__)*"/g163303.03170_CHEASE"
const free_geqdsk = (@__DIR__)*"/g163303.03170_EFIT"
const good_currents = [-197968.56551716285, -175313.42410039817, -52907.88075466734, -256837.87396872052, 485346.0233782668, -214060.41838723494, -151469.54405769982, -330872.8548105091, -488.24739294769034, -194334.88010685574, -136597.64936162502, -159445.93012418432, 16601.967623453867, 53881.80298066415, -237456.52831624495, -230790.53129979773, 85387.28458689603, -23680.644329152925]

@testset "coil_current" begin
    g0 = readg(fixed_geqdsk)
    cc0 = cocos(g0,clockwise_phi=false).cocos

    ccs = collect([1:8;11:18])
    ccp = rand(ccs)
    for cc in ccs
        gfixed = transform_cocos(g0, cc0, cc)
        Gfixed = efit(gfixed, cc)
    
        currents = fixed_eq_currents(Gfixed,coils_D3D)
        @test currents ≈ Gfixed.cocos.sigma_RpZ * good_currents

        if cc == ccp
            println("Plotting for COCOS ", cc)
            gfree = transform_cocos(readg(free_geqdsk), cc0, cc)
            Gfree = efit(gfree, cc)
            check_fixed_eq_currents(Gfixed,coils_D3D,currents,Gfree)
        end
    end
end

@testset "Solovev_current" begin
    δ = 0.5
    ϵ = 0.32
    κ = 1.7
    B0 = 2.0
    R0 = 1.8
    qstar = 1.57
    alpha = -0.155
    S = solovev(B0, R0, ϵ, δ, κ, alpha, qstar,B0_dir=1,Ip_dir=1)
    currents = fixed_eq_currents(S,coils_D3D)
    check_fixed_eq_currents(S,coils_D3D,currents)
end