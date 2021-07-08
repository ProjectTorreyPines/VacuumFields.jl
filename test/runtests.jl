using AD_GS
using Test
using EFIT
using Equilibrium
#using DelimitedFiles

const active = [
    "current_cocos",
    "current_BtIp",
    "current_Solovev"
]

const coils = coils_D3D

if "current_cocos" in active
    @testset "current_cocos" begin
        fixed_geqdsk = (@__DIR__)*"/equilibria/g163303.03170_fix"
        free_geqdsk  = (@__DIR__)*"/equilibria/g163303.03170_free"

        if coils == coils_D3D_points
            good_currents = [-197968.56551716285, -175313.42410039817, -52907.88075466734,
                            -256837.87396872052, 485346.0233782668, -214060.41838723494,
                            -151469.54405769982, -330872.8548105091, -488.24739294769034,
                            -194334.88010685574, -136597.64936162502, -159445.93012418432,
                            16601.967623453867, 53881.80298066415, -237456.52831624495,
                            -230790.53129979773, 85387.28458689603, -23680.644329152925]
        elseif coils == coils_D3D
            good_currents = [-196476.6455553291, -180470.6383359014, -54474.37322174671,
                             -201354.39635394656, 364183.31928296294, -213850.00567157066,
                             -152910.2685568697, -277343.37135676463, 2149.816576966399,
                             -199728.91752482625, -129092.5550445522, -167017.5580650013,
                             -1401.6721399478631, 117437.4109138634, -237058.76256911657,
                             -232926.81215186242, 47000.47824298689, -16537.43992243954]
        end

        g0 = readg(fixed_geqdsk)
        cc0 = cocos(g0,clockwise_phi=false).cocos

        ccs = collect([1:8;11:18])
        ccp = rand(ccs)
        for cc in ccs
            gfixed = transform_cocos(g0, cc0, cc)
            Gfixed = efit(gfixed, cc)

            currents = fixed_eq_currents(Gfixed,coils)
            @test currents ≈ Gfixed.cocos.sigma_RpZ * good_currents

            if cc == ccp
                println("Plotting for COCOS ", cc)
                gfree = transform_cocos(readg(free_geqdsk), cc0, cc)
                Gfree = efit(gfree, cc)
                check_fixed_eq_currents(Gfixed,coils,currents,Gfree)
            end
        end
    end
end

if "current_BtIp" in active
    @testset "current_BtIp" begin
        fixed_pp = (@__DIR__)*"/equilibria/g150219.03200_fix"
        free_pp  = (@__DIR__)*"/equilibria/g150219.03200_free"
    
        fixed_pm = (@__DIR__)*"/equilibria/g159177.02700_fix"
        free_pm  = (@__DIR__)*"/equilibria/g159177.02700_free"

        fixed_mp = (@__DIR__)*"/equilibria/g133221.01151_fix"
        free_mp  = (@__DIR__)*"/equilibria/g133221.01151_free"

        fixed_mm = (@__DIR__)*"/equilibria/g153298.04400_fix"
        free_mm  = (@__DIR__)*"/equilibria/g153298.04400_free"

        #C = readdlm((@__DIR__)*"/equilibria/currents.txt")

        EQs_fixed = [fixed_pp,fixed_pm,fixed_mp,fixed_mm]
        EQs_free  = [free_pp, free_pm, free_mp, free_mm]

        for i in 1:4
            gfixed = readg(EQs_fixed[i])
            cc = cocos(gfixed, clockwise_phi=false).cocos
            Gfixed = efit(gfixed, cc)
            currents = fixed_eq_currents(Gfixed,coils)
            # p = plot(C[i,:])
            # plot!(currents)
            # display(p)
            #@test currents ≈ Gfixed.cocos.sigma_RpZ * good_currents
            gfree = readg(EQs_free[i])
            Gfree = efit(gfree, cc)
            check_fixed_eq_currents(Gfixed,coils,currents,Gfree)
        end
    end
end

if "current_Solovev" in active
    @testset "current_Solovev" begin
        δ = 0.5
        ϵ = 0.32
        κ = 1.7
        B0 = 2.0
        R0 = 1.8
        qstar = 1.57
        alpha = -0.155
        S = solovev(B0, R0, ϵ, δ, κ, alpha, qstar,B0_dir=1,Ip_dir=1)
        currents = fixed_eq_currents(S,coils)
        check_fixed_eq_currents(S,coils,currents)
    end
end