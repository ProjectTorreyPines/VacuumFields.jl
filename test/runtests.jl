using AD_GS
using Test
using EFIT
using Equilibrium
import DelimitedFiles
using Plots

const active = [
    #"current_cocos",
    #"current_BtIp",
    #"current_Solovev",
    "current_breakdown"
]

# Change to coils_D3D_points to run with singular coils
const coils = coils_D3D_points

if "current_cocos" in active
    @testset "current_cocos" begin
        fixed_geqdsk = (@__DIR__)*"/equilibria/g163303.03170_fix"
        free_geqdsk  = (@__DIR__)*"/equilibria/g163303.03170_free"

        # These are the currents computed by EFIT, without the E-coils included
        eqdsk_currents = [-118002.985, -95272.9804, 9404.27229, 
                           45974.8807, 27840.9167, -206295.043,
                          -141605.015, 13459.9473, -1087.47886,
                          -113719.456, -59571.7105, -69097.6062,
                           98250.1947, 143905.727, -230652.605,
                           -225334.449, 116828.799, 8614.80714]

        if coils == coils_D3D_points
            good_currents = [-116235.28907610463, -98586.2314474624, 65689.70513193491,
                             -469345.0846858057, 1.429593664421151e6, -204968.27495366192,
                             -145109.5456392649, -868933.9278994281, 119612.3742593849,
                             -117495.15065519481, -56048.70738466947, -78894.1800648648,
                              108238.61937784836, 159156.2152710073, -230654.60093398875,
                             -221695.36575228465, 89848.68347638033, 13032.984674759902]
        elseif coils == coils_D3D
            good_currents = [-113885.7197764252, -105436.95737589894, 60341.419031375575,
                             -340056.74903837184, 1.1179306735975991e6, -204783.63904252168,
                             -146579.64500224337, -711735.8837820432, 115364.7725883161,
                             -123685.36081361726, -47752.21100839885, -89037.91895772259,
                              98446.71882299529, 211077.24485467927, -230247.60310567988,
                             -223867.46897678653, 58818.18486954898, 18730.207514097914]
        end

        gfixed = readg(fixed_geqdsk)
        cc0 = cocos(gfixed,clockwise_phi=false).cocos
        gfree = readg(free_geqdsk)

        ccs = collect([1:8;11:18])
        ccp = rand(ccs)

        for cc in ccs
            Gfixed = efit(transform_cocos(gfixed, cc0, cc), cc)
            Gfree = efit(transform_cocos(gfree, cc0, cc), cc)
            _,ψbound = psi_limits(Gfree)
            currents = fixed_eq_currents(Gfixed,coils,ψbound)
            @test currents ≈ Gfixed.cocos.sigma_RpZ * good_currents

            if cc == cc0
                p = plot(eqdsk_currents, linewidth=3, linecolor=:black)
                plot!(currents, linewidth=3, linecolor=:blue)
                display(p)
            end

            if cc == ccp
                println("Plotting for COCOS ", cc)
                p = check_fixed_eq_currents(Gfixed,coils,currents,Gfree)
                display(p)
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

        C = DelimitedFiles.readdlm((@__DIR__)*"/equilibria/currents.txt")

        EQs_fixed = [fixed_pp,fixed_pm,fixed_mp,fixed_mm]
        EQs_free  = [free_pp, free_pm, free_mp, free_mm]

        for i in 1:4
            gfixed = readg(EQs_fixed[i])
            cc = cocos(gfixed, clockwise_phi=false).cocos
            Gfixed = efit(gfixed, cc)

            gfree = readg(EQs_free[i])
            Gfree = efit(gfree, cc)
            _,ψbound = psi_limits(Gfree)

            #  Currents from EFIT
            p = plot(C[i,:], linewidth=3, linecolor=:black)

            # Currents with ψbound=0
            c0 = fixed_eq_currents(Gfixed, coils, 0.0)
            plot!(c0, linewidth=3, linecolor=:red)

            # Currents with ψbound from EFIT
            cb = fixed_eq_currents(Gfixed, coils, ψbound)
            plot!(cb, linewidth=3, linecolor=:blue)

            # Currents minimized
            # cm = fixed_eq_currents(Gfixed, coils, ψbound, minimize_currents=0.01)
            # plot!(currents)
            display(p)

            # p = check_fixed_eq_currents(Gfixed,coils,c0,Gfree)
            # display(p)
            # p = check_fixed_eq_currents(Gfixed,coils,cb,Gfree)
            # display(p)

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
        p = check_fixed_eq_currents(S,coils,currents)
        display(p)
    end
end

if "current_breakdown" in active
    @testset "current_breakdown" begin

        # OH coils
        Noh = 100
        Roh = 0.8*ones(Noh)
        Zoh = range(-1.0, 1.0, length=Noh)
        coils_OH = collect(zip(Roh,Zoh))
        currents_OH = 1000.0*ones(Noh)

        # PF coils
        # From DIII-D, but remove coils along CS
        coils_PF = [(1.0041, 1.5169), (2.6124, 0.4376),
                    (2.3733, 1.1171), (1.2518, 1.6019),
                    (1.689, 1.5874), (1.0025, -1.5169),
                    (2.6124, -0.4376), (2.3834, -1.1171),
                    (1.2524, -1.6027), (1.6889, -1.578)]

        # boundary of field null
        # arbitrary ellipse
        θ = range(-π,π,length=129)
        R₀ = 1.4
        Z₀ = 0.0
        a₀ = 0.4
        ϵ = 1.0
        Rp = R₀ .+ a₀.*cos.(θ)
        Zp = Z₀ .+ ϵ*a₀.*sin.(θ)

        # ψ we want on this boundary
        # HOW DOES THIS GET DETERMINED?
        # Probably the highest \psi from the OH coils?
        ψp = -0.015

        # Find PF coil currents to make field null
        currents_PF = currents_to_match_ψp(1.0, ψp, Rp, Zp, coils_PF,
                                           fixed_coils=coils_OH,
                                           fixed_currents=currents_OH)

        println(1e-3*currents_PF)

        # Plot flux from all coils
        coils_all = [coils_PF; coils_OH]
        currents_all = [currents_PF; currents_OH]
        p = plot_coil_flux(1.0,coils_all,currents_all,ψp,clim=(-0.03,-0.00),resolution=1025,Rmin=0.5,Rmax=3.0,Zmin=-1.5,Zmax=1.5)
        plot!(p,Rp,Zp)
        display(p)
    end
end
