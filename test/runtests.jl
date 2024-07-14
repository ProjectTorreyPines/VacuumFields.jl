using VacuumFields
using Test
using EFIT
using MXHEquilibrium
import DelimitedFiles
using Plots

const active = [
    "current_cocos",
    "current_BtIp",
    "current_Solovev",
    "current_breakdown",
    "current_xpoint"
]

const coils_D3D_matrix = [
    [8.6080e-01 1.6830e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6140e-01 5.0810e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6280e-01 8.4910e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6110e-01 1.1899e+00 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [1.0041e+00 1.5169e+00 1.3920e-01 1.1940e-01 4.5000e+01 9.0000e+01]
    [2.6124e+00 4.3760e-01 1.7320e-01 1.9460e-01 0.0000e+00 9.2400e+01]
    [2.3733e+00 1.1171e+00 1.8800e-01 1.6920e-01 0.0000e+00 1.0806e+02]
    [1.2518e+00 1.6019e+00 2.3490e-01 8.5100e-02 0.0000e+00 9.0000e+01]
    [1.6890e+00 1.5874e+00 1.6940e-01 1.3310e-01 0.0000e+00 9.0000e+01]
    [8.6080e-01 -1.7370e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6070e-01 -5.1350e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6110e-01 -8.5430e-01 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [8.6300e-01 -1.1957e+00 5.0800e-02 3.2110e-01 0.0000e+00 9.0000e+01]
    [1.0025e+00 -1.5169e+00 1.3920e-01 1.1940e-01 -4.5000e+01 9.0000e+01]
    [2.6124e+00 -4.3760e-01 1.7320e-01 1.9460e-01 0.0000e+00 -9.2400e+01]
    [2.3834e+00 -1.1171e+00 1.8800e-01 1.6920e-01 0.0000e+00 -1.0806e+02]
    [1.2524e+00 -1.6027e+00 2.3490e-01 8.5100e-02 0.0000e+00 9.0000e+01]
    [1.6889e+00 -1.5780e+00 1.6940e-01 1.3310e-01 0.0000e+00 9.0000e+01]]

const coils_D3D_points = [PointCoil(coils_D3D_matrix[i, 1], coils_D3D_matrix[i, 2]) for i in 1:size(coils_D3D_matrix)[1]]
const coils_D3D = [ParallelogramCoil(coils_D3D_matrix[i, :]...) for i in 1:size(coils_D3D_matrix)[1]]

# Change to coils_D3D_points to run with singular coils
const coils = coils_D3D

if "current_cocos" in active
    @testset "current_cocos" begin
        fixed_geqdsk = (@__DIR__) * "/equilibria/g163303.03170_fix"
        free_geqdsk = (@__DIR__) * "/equilibria/g163303.03170_free"

        # These are the currents computed by EFIT, without the E-coils included
        eqdsk_currents = [
            -118002.985, -95272.9804, 9404.27229,
            45974.8807, 27840.9167, -206295.043,
            -141605.015, 13459.9473, -1087.47886,
            -113719.456, -59571.7105, -69097.6062,
            98250.1947, 143905.727, -230652.605,
            -225334.449, 116828.799, 8614.80714]

        if coils == coils_D3D_points
            good_currents = [
                -116235.28907610463, -98586.2314474624, 65689.70513193491,
                -469345.0846858057, 1.429593664421151e6, -204968.27495366192,
                -145109.5456392649, -868933.9278994281, 119612.3742593849,
                -117495.15065519481, -56048.70738466947, -78894.1800648648,
                108238.61937784836, 159156.2152710073, -230654.60093398875,
                -221695.36575228465, 89848.68347638033, 13032.984674759902]
        elseif coils == coils_D3D
            good_currents_sd = [
                -123635.83053274715, -101040.44981352903, -8692.175095363986,
                99456.30564584397, 44148.55159103844, -212076.3721057026,
                -136048.95577378047, -41156.17734011926, 14409.220391036204,
                -131177.93090524126, -47934.92152827256, -93944.56599139964,
                119624.81517763663, 164320.37132246612, -235830.7658811443,
                -221515.04511039687, 97744.20420786832, 11170.040045183569]
            good_currents_dd = [
                -123796.43108312064, -101086.71033357759, 6330.115901381709,
                -55169.25566862989, 512296.8855078304, -211651.50879072546,
                -139010.09620696527, -359186.05410263385, 65945.91628320269,
                -131121.50001405698, -48038.765505325864, -92503.48058979367,
                110876.6061240047, 186605.5564113924, -235893.39824374302,
                -221782.0227385224, 82690.44538168146, 14080.376817403361]
        end

        gfixed = readg(fixed_geqdsk)
        cc0 = cocos(gfixed; clockwise_phi=false).cocos
        gfree = readg(free_geqdsk)

        ccs = collect([1:8; 11:18])
        ccp = rand(ccs)

        for cc in ccs
            println("Testing for COCOS ", cc)
            Gfixed = efit(transform_cocos(gfixed, cc0, cc), cc)
            Gfree = efit(transform_cocos(gfree, cc0, cc), cc)
            _, ψbound = psi_limits(Gfree)
            currents = fixed_eq_currents(Gfixed, coils, AbstractCoil[], ψbound)
            if cc < 10
                good_currents = good_currents_sd
            else
                good_currents = good_currents_dd
            end

            if cc == cc0
                p = plot(eqdsk_currents; linewidth=3, linecolor=:black)
                plot!(currents; linewidth=3, linecolor=:blue)
                display(p)
            end

            if cc == ccp
                println("Plotting for COCOS ", cc)
                p = check_fixed_eq_currents(Gfixed, coils, Gfree)
                display(p)
            end

            @test currents ≈ Gfixed.cocos.sigma_RpZ * good_currents
        end
    end
end

if "current_BtIp" in active
    @testset "current_BtIp" begin
        fixed_pp = (@__DIR__) * "/equilibria/g150219.03200_fix"
        free_pp = (@__DIR__) * "/equilibria/g150219.03200_free"

        fixed_pm = (@__DIR__) * "/equilibria/g159177.02700_fix"
        free_pm = (@__DIR__) * "/equilibria/g159177.02700_free"

        fixed_mp = (@__DIR__) * "/equilibria/g133221.01151_fix"
        free_mp = (@__DIR__) * "/equilibria/g133221.01151_free"

        fixed_mm = (@__DIR__) * "/equilibria/g153298.04400_fix"
        free_mm = (@__DIR__) * "/equilibria/g153298.04400_free"

        C = DelimitedFiles.readdlm((@__DIR__) * "/equilibria/currents.txt")

        EQs_fixed = [fixed_pp, fixed_pm, fixed_mp, fixed_mm]
        EQs_free = [free_pp, free_pm, free_mp, free_mm]

        for i in 1:4
            gfixed = readg(EQs_fixed[i])
            cc = cocos(gfixed; clockwise_phi=false).cocos
            Gfixed = efit(gfixed, cc)

            gfree = readg(EQs_free[i])
            Gfree = efit(gfree, cc)
            _, ψbound = psi_limits(Gfree)

            #  Currents from EFIT
            p = plot(C[i, :]; linewidth=3, linecolor=:black)

            # Currents with ψbound=0
            c0 = fixed_eq_currents(Gfixed, coils, AbstractCoil[], 0.0)
            plot!(c0; linewidth=3, linecolor=:red)

            # Currents with ψbound from EFIT
            cb = fixed_eq_currents(Gfixed, coils, AbstractCoil[], ψbound)
            plot!(cb; linewidth=3, linecolor=:blue)

            # Currents minimized
            # cm = fixed_eq_currents(Gfixed, coils, AbstractCoil[], ψbound, minimize_currents=0.01)
            # plot!(currents)
            display(p)

            # p = check_fixed_eq_currents(Gfixed,coils,Gfree)
            # display(p)
            # p = check_fixed_eq_currents(Gfixed,coils,Gfree)
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
        S = solovev(B0, R0, ϵ, δ, κ, alpha, qstar; B0_dir=1, Ip_dir=1)
        currents = fixed_eq_currents(S, coils)
        p = check_fixed_eq_currents(S, coils)
        display(p)
    end
end

if "current_breakdown" in active
    @testset "current_breakdown" begin

        # OH coils
        coils_OH = coils_D3D[[1, 2, 3, 4, 10, 11, 12, 13]]
        for coil in coils_OH
            coil.current = 10000.0
        end

        # PF coils
        # From DIII-D, but remove coils along CS
        coils_PF = coils_D3D[[5, 6, 7, 8, 9, 14, 15, 16, 17, 18]]

        # boundary of field null
        # arbitrary ellipse
        θ = range(-π, π; length=129)
        R₀ = 1.6
        Z₀ = 0.0
        a₀ = 0.3
        ϵ = 1.5
        Rp = R₀ .+ a₀ .* cos.(θ)
        Zp = Z₀ .+ ϵ * a₀ .* sin.(θ)

        # ψ we want on this boundary
        # HOW DOES THIS GET DETERMINED?
        # Probably the highest \psi from the OH coils?
        ψp_constant = -0.012

        # account for effect of fixed coils on ψp_constant
        Bp_fac, ψp, Rp, Zp = field_null_on_boundary(ψp_constant, Rp, Zp, coils_OH)

        # Find PF coil currents to make field null
        currents_PF = currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils_PF)

        # Plot flux from all coils
        coils_all = vcat(coils_PF, coils_OH)
        p = plot_coil_flux(Bp_fac, coils_all, ψp_constant; clim=(2.0 * ψp_constant, 0.0), Rmin=0.5, Rmax=3.0, Zmin=-1.8, Zmax=1.8)
        for coil in coils_all
            plot_coil(coil)
        end
        plot!(p, Rp, Zp)
        display(p)
    end
end

if "current_xpoint" in active
    @testset "current_xpoint" begin

        fixed_geqdsk = (@__DIR__) * "/equilibria/g163303.03170_fix"
        free_geqdsk = (@__DIR__) * "/equilibria/g163303.03170_free"
        gfixed = readg(fixed_geqdsk)
        gfree = readg(free_geqdsk)
        cc = cocos(gfixed; clockwise_phi=false).cocos
        Gfixed = efit(gfixed, cc)
        Gfree = efit(gfree, cc)
        _, ψbound = psi_limits(Gfree)

        Rx = [1.7]
        Zx = [-1.5]
        Bp_fac, ψp, Rp, Zp = ψp_on_fixed_eq_boundary(Gfixed, AbstractCoil[], ψbound; Rx=Rx, Zx=Zx)

        weights = ones(length(Rp))
        weights[end] = 1.0

        # X points can't get added later
        # Outside the LCFS, there's still a flux contribution
        # from the image currents (IS THAT TRUE???)

        currents = currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; weights=weights, λ_regularize=1E-14)

        p = check_fixed_eq_currents(Gfixed, coils, Gfree; Rmin=0.5, Rmax=3.0, Zmin=-2.0, Zmax=2.0)
        display(p)
    end
end
