using VacuumFields
using Test
using MXHEquilibrium
import EFIT: readg
import DelimitedFiles
using Plots
using MAT

const Ip = 1e5
const resistance = 1e-6
const turns = 2

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

const coils = [QuadCoil(ParallelogramCoil(pg..., Ip; resistance, turns)) for pg in eachrow(coils_D3D_matrix)]


function load_par(params)
    R = params[2]
    Z = params[1]
    ΔR = params[4]
    ΔZ = params[3]
    θ₁ = params[5]
    #θ₂ = params[6] - 90.0
    # ??? What's the equation here?
    θ₂ = (params[6] == 0.0) ? 90.0 : params[6]
    return R, Z, ΔR, ΔZ, θ₁, θ₂
end

@testset "VacuumFields.jl" begin

    @test PointCoil(coils_D3D_matrix[1, 1], coils_D3D_matrix[1, 2], Ip; resistance, turns) isa PointCoil
    parcoil = ParallelogramCoil(coils_D3D_matrix[1, :]..., Ip; resistance, turns)
    @test parcoil isa ParallelogramCoil
    @test QuadCoil(parcoil) isa QuadCoil
    qcoil = QuadCoil([1.0, 2.0, 3.0, 1.0], [-1.0, -1.0, 2.0, 1.0], Ip; resistance, turns)
    @test VacuumFields.area(qcoil) ≈ 3.5
    R, Z = 10.0, 10.0
    @test 2π * VacuumFields.μ₀ * VacuumFields.Green(qcoil, R, Z) * VacuumFields.VacuumFields.current(qcoil) ≈ VacuumFields.ψ(qcoil, R, Z)
    @test DistributedCoil(parcoil) isa DistributedCoil

    mcoil = MultiCoil(coils)
    Nc = length(coils)
    @test VacuumFields.ψ(mcoil, 2.0, 0.0) ≈ sum(VacuumFields.ψ(coil, 2.0, 0.0) for coil in coils)
    @test VacuumFields.current(mcoil) ≈ Nc * Ip
    @test VacuumFields.resistance(mcoil) ≈ Nc * resistance
    @test VacuumFields.turns(mcoil) == Nc * turns
    Ic1 = VacuumFields.current(coils[1])
    VacuumFields.set_current!(mcoil, Ip)
    @test VacuumFields.current(mcoil.coils[1]) ≈ Ip / Nc
    VacuumFields.set_current!(mcoil, Nc * Ip)
    @test VacuumFields.current(mcoil.coils[1]) ≈ Ip

    geqdsk = readg((@__DIR__) * "/equilibria/g184250.01740"; set_time=0.0)
    cc0 = cocos(geqdsk; clockwise_phi=false).cocos
    EQ = efit(transform_cocos(geqdsk, cc0, 11), 11)

    vars = matread((@__DIR__) * "/equilibria/d3d_vs_eq3_.mat")
    matcoils = [ParallelogramCoil(load_par(vars["fcdata"][:,k])..., vars["Ic"][2+k] * Int(vars["fcnturn"][k]); resistance=vars["Rc"][k], turns=Int(vars["fcnturn"][k])) for k in eachindex(vars["fcnturn"])]
    append!(matcoils, ParallelogramCoil(load_par(vars["vvdata"][:,k])...; resistance=vars["Rv"][k]) for k in eachindex(vars["Rv"]))

    @test isapprox(stability_margin(EQ, matcoils), vars["stability"]["m_s_noe"]; rtol=5e-2)
    γ, τ, _ = normalized_growth_rate(EQ, matcoils)
    γ_ts, τ_ts = vars["stability"]["gamma_massless_noe"], vars["stability"]["tauv_massless"]
    @test isapprox(γ, γ_ts; rtol=5e-2)
    @test isapprox(τ, τ_ts; rtol=5e-2)
end



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


    good_currents_sd = [
        -125255.22802038345, -85732.61193982523, 6908.577576682961,
        29275.365128753736, 25359.985822088915, -211204.251865932,
        -134448.83811184048, 23488.525199361036, -5294.195981690289,
        -107043.98891397915, -75142.7511909334, -48739.08096675255,
        77599.0070426326, 150830.88444440905, -221402.16070576818,
        -228650.7045731426, 123465.62290190649, 926.6664242527331
        ]
    good_currents_dd = [
        -115274.78058122244, -100066.62138053821, 16066.723028895212,
        31472.34544931911, 30233.29453776579, -212811.09342069676,
        -131915.01608336507, 25955.390168425336, -12890.287314473018,
        -116340.44880288554, -64577.55852943135, -47485.319472859,
        31959.71428826521, 249319.4096091399, -217067.7110646716,
        -240216.66971257926, 46785.55843708216, 33063.62059136332
        ]

    gfixed = readg(fixed_geqdsk; set_time=0.0)
    cc0 = cocos(gfixed; clockwise_phi=false).cocos
    gfree = readg(free_geqdsk; set_time=0.0)

    ccs = collect([1:8; 11:18])

    ix = argmin(gfree.zbbbs)
    Rx, Zx = gfree.rbbbs[ix], gfree.zbbbs[ix]

    for cc in ccs
        println("Testing for COCOS ", cc)
        Gfixed = efit(transform_cocos(gfixed, cc0, cc), cc)
        Gfree = efit(transform_cocos(gfree, cc0, cc), cc)
        _, ψbound = psi_limits(Gfree)
        flux_cps = boundary_control_points(Gfree, 0.999, ψbound)
        saddle_cps = [SaddleControlPoint(Rx, Zx)]
        currents, _ = find_coil_currents!(coils, Gfixed; flux_cps, saddle_cps, λ_regularize=1e-14)
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

        @test currents ≈ Gfixed.cocos.sigma_RpZ * good_currents
    end
end

@testset "fixed2free" begin
    fixed_geqdsk = (@__DIR__) * "/equilibria/g163303.03170_fix"
    free_geqdsk = (@__DIR__) * "/equilibria/g163303.03170_free"

    gfixed = readg(fixed_geqdsk; set_time=0.0)
    cc0 = cocos(gfixed; clockwise_phi=false).cocos
    gfree = readg(free_geqdsk; set_time=0.0)
    Gfixed = efit(transform_cocos(gfixed, cc0, 11), 11)
    Gfree = efit(transform_cocos(gfree, cc0, 11), 11)

    _, ψbound = psi_limits(Gfree)
    Rs, Zs = gfree.r, gfree.z

    ix = argmin(gfree.zbbbs)
    Rx, Zx = gfree.rbbbs[ix], gfree.zbbbs[ix]

    flux_cps = boundary_control_points(Gfree, 0.999, ψbound)
    saddle_cps = [SaddleControlPoint(Rx, Zx)]
    fixed2free(Gfixed, coils, Rs, Zs; flux_cps, saddle_cps)
end

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
        gfixed = MXHEquilibrium.readg(EQs_fixed[i]; set_time=0.0)
        cc = MXHEquilibrium.cocos(gfixed; clockwise_phi=false).cocos
        Gfixed = MXHEquilibrium.efit(gfixed, cc)

        gfree = MXHEquilibrium.readg(EQs_free[i]; set_time=0.0)
        Gfree = MXHEquilibrium.efit(gfree, cc)
        _, ψbound = psi_limits(Gfree)

        #  Currents from EFIT
        p = plot(C[i, :]; linewidth=3, linecolor=:black)

        # Currents with ψbound=0
        flux_cps = boundary_control_points(Gfixed, 0.999, 0.0)
        c0, _ = find_coil_currents!(coils, Gfixed; flux_cps, λ_regularize=1e-14)

        plot!(c0; linewidth=3, linecolor=:red)

        # Currents with ψbound from EFIT
        flux_cps = boundary_control_points(Gfixed, 0.999, ψbound)
        cb, _ = find_coil_currents!(coils, Gfixed; flux_cps, λ_regularize=1e-14)
        plot!(cb; linewidth=3, linecolor=:blue)
        display(p)

    end
end

# @testset "current_breakdown" begin

#     # OH coils
#     coils_OH = coils_D3D[[1, 2, 3, 4, 10, 11, 12, 13]]
#     for coil in coils_OH
#         coil.current = 10000.0
#     end

#     # PF coils
#     # From DIII-D, but remove coils along CS
#     coils_PF = coils_D3D[[5, 6, 7, 8, 9, 14, 15, 16, 17, 18]]

#     # boundary of field null
#     # arbitrary ellipse
#     θ = range(-π, π; length=129)
#     R₀ = 1.6
#     Z₀ = 0.0
#     a₀ = 0.3
#     ϵ = 1.5
#     Rp = R₀ .+ a₀ .* cos.(θ)
#     Zp = Z₀ .+ ϵ * a₀ .* sin.(θ)

#     # ψ we want on this boundary
#     # HOW DOES THIS GET DETERMINED?
#     # Probably the highest \psi from the OH coils?
#     ψp_constant = -0.012

#     # account for effect of fixed coils on ψp_constant
#     Bp_fac, ψp, Rp, Zp = field_null_on_boundary(ψp_constant, Rp, Zp, coils_OH)

#     # Find PF coil currents to make field null
#     currents_PF = currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils_PF)

#     # Plot flux from all coils
#     coils_all = vcat(coils_PF, coils_OH)
#     p = plot_coil_flux(Bp_fac, coils_all, ψp_constant; clim=(2.0 * ψp_constant, 0.0), Rmin=0.5, Rmax=3.0, Zmin=-1.8, Zmax=1.8)
#     for coil in coils_all
#         plot_coil(coil)
#     end
#     plot!(p, Rp, Zp)
#     display(p)
# end
