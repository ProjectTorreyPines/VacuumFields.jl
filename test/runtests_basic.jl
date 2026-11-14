@testset "VacuumFields.jl" begin

    @test PointCoil(coils_D3D_matrix[1, 1], coils_D3D_matrix[1, 2], Icpt; resistance, turns) isa PointCoil
    parcoil = ParallelogramCoil(coils_D3D_matrix[1, :]..., Icpt; resistance, turns)
    @test parcoil isa ParallelogramCoil
    @test QuadCoil(parcoil) isa QuadCoil
    qcoil = QuadCoil([1.0, 2.0, 3.0, 1.0], [-1.0, -1.0, 2.0, 1.0], Icpt; resistance, turns)
    @test VacuumFields.area(qcoil) ≈ 3.5
    R, Z = 10.0, 10.0
    @test 2π * VacuumFields.μ₀ * VacuumFields.Green(qcoil, R, Z) * VacuumFields.VacuumFields.current_per_turn(qcoil) ≈ VacuumFields.ψ(qcoil, R, Z)
    @test DistributedCoil(parcoil) isa DistributedCoil

    orientation = ones(Int, length(coils))
    orientation[1:2:end] .= -1
    mcoil = MultiCoil(coils, orientation)
    VacuumFields.set_current_per_turn!(mcoil, Icpt)
    @test VacuumFields.current_per_turn(mcoil) ≈ Icpt
    Nc = length(coils)
    r, z = 2.0, 0.0
    @test VacuumFields.ψ(mcoil, r, z) ≈ sum(VacuumFields.ψ(coil, r, z) * mcoil.orientation[k] for (k, coil) in enumerate(coils)) # these all have Icpt current/turn
    @test VacuumFields.ψ(mcoil, r, z) ≈ sum(VacuumFields.ψ(coil, r, z) for (k, coil) in enumerate(mcoil.coils)) # these have Icpt current/turn with appropriate sign
    @test VacuumFields.ψ(mcoil, r, z) ≈ 2π * VacuumFields.μ₀ * Icpt * sum(VacuumFields.Green(coil, r, z) * mcoil.orientation[k] for (k, coil) in enumerate(mcoil.coils))
    @test VacuumFields.ψ(mcoil, r, z) ≈ 2π * VacuumFields.μ₀ * Icpt * VacuumFields.Green(mcoil, r, z)
    @test VacuumFields.resistance(mcoil) ≈ Nc * resistance
    @test VacuumFields.turns(mcoil) == Nc * turns
    VacuumFields.set_current_per_turn!(mcoil, 2.0 * Icpt)
    @test VacuumFields.current_per_turn(mcoil.coils[1]) * mcoil.orientation[1] ≈ 2.0 * Icpt

    geqdsk = readg((@__DIR__) * "/equilibria/g184250.01740"; set_time=0.0)
    cc0 = cocos(geqdsk; clockwise_phi=false).cocos
    EQ = efit(transform_cocos(geqdsk, cc0, 11), 11)

    vars = matread((@__DIR__) * "/equilibria/d3d_vs_eq3_.mat")
    matcoils = [ParallelogramCoil(load_par(vars["fcdata"][:,k])..., vars["Ic"][2+k]; resistance=vars["Rc"][k], turns=Int(vars["fcnturn"][k])) for k in eachindex(vars["fcnturn"])]
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


    gfixed = readg(fixed_geqdsk; set_time=0.0)
    cc0 = cocos(gfixed; clockwise_phi=false).cocos
    gfree = readg(free_geqdsk; set_time=0.0)

    ccs = collect([1:8; 11:18])

    ix = argmin(gfree.zbbbs)
    Rx, Zx = gfree.rbbbs[ix], gfree.zbbbs[ix]
    saddle_cps = [SaddleControlPoint(Rx, Zx)]
    for cc in ccs
        println("Testing for COCOS ", cc)
        Gfixed = efit(transform_cocos(gfixed, cc0, cc), cc)
        Gfree = efit(transform_cocos(gfree, cc0, cc), cc)
        _, ψbound = psi_limits(Gfree)
        flux_cps = boundary_control_points(Gfree, 0.999, ψbound)
        currents_per_turn, cost = find_coil_currents!(coils, Gfixed; flux_cps, saddle_cps, λ_regularize=1e-14)
        currents = currents_per_turn .* VacuumFields.turns.(coils)
        if cc == cc0
            p = plot(eqdsk_currents; linewidth=3, linecolor=:black)
            plot!(currents; linewidth=3, linecolor=:blue)
            display(p)
        end

        # The fit is nearly degenerate in the currents (tiny regularization), so the individual currents
        # depend on round-off in the linear solve. Test what the fit is for instead:
        # the control points are matched ...
        @test cost < 1e-2
        # ... by currents that follow those of EFIT, with the sign of the COCOS
        ref = Gfixed.cocos.sigma_RpZ * eqdsk_currents
        @test sum(currents .* ref) / sqrt(sum(abs2, currents) * sum(abs2, ref)) > 0.9
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
        cpt0, _ = find_coil_currents!(coils, Gfixed; flux_cps, λ_regularize=1e-14)
        c0 = cpt0 .* VacuumFields.turns.(coils)

        plot!(c0; linewidth=3, linecolor=:red)

        # Currents with ψbound from EFIT
        flux_cps = boundary_control_points(Gfixed, 0.999, ψbound)
        cptb, _ = find_coil_currents!(coils, Gfixed; flux_cps, λ_regularize=1e-14)
        cb = cptb .* VacuumFields.turns.(coils)
        plot!(cb; linewidth=3, linecolor=:blue)
        display(p)

    end
end

@testset "all cp types" begin
    tol = 100 * sqrt(eps())
    θ = range(0, 2π, 19)[1:end-1]
    Rs = 1.5 .+ 0.25 * cos.(θ)
    Zs = 0.5 * sin.(θ)

    # 18 control points total
    flux_cps = [FluxControlPoint(Rs[3k+1], Zs[3k+1], 1.0 * k) for k in 0:3] # 4 flux control points
    saddle_cps = [SaddleControlPoint(Rs[3k+2], Zs[3k+2]) for k in 0:2:4] # 3 saddle control points (2 constaints each)
    iso_cps = [IsoControlPoint(Rs[3k+3], Zs[3k+3], Rs[1], Rs[2], k - 1.0) for k in 0:3] # 4 iso control points with offsets
    field_cps = [FieldControlPoint(Rs[3k+1], Zs[3k+1], 1.0 * k, 1.0 * k * k) for k in 0:3] # 4 field control points
    currents_per_turn, _ = find_coil_currents!(coils, nothing; flux_cps, saddle_cps, iso_cps, field_cps)
    for flux_cp in flux_cps
        psi = sum(VacuumFields.ψ(coil, flux_cp.R, flux_cp.Z) for coil in coils)
        @test isapprox(psi, flux_cp.target; atol=tol)
    end
    for saddle_cp in saddle_cps
        dpsidr = sum(VacuumFields.dψ_dR(coil, saddle_cp.R, saddle_cp.Z) for coil in coils)
        @test isapprox(dpsidr, 0.0; atol=tol)
        dpsidz = sum(VacuumFields.dψ_dZ(coil, saddle_cp.R, saddle_cp.Z) for coil in coils)
        @test isapprox(dpsidz, 0.0; atol=tol)
    end
    for iso_cp in iso_cps
        psi1 = sum(VacuumFields.ψ(coil, iso_cp.R1, iso_cp.Z1) for coil in coils)
        psi2 = sum(VacuumFields.ψ(coil, iso_cp.R2, iso_cp.Z2) for coil in coils)
        @test isapprox(psi1 - psi2, iso_cp.offset; atol=tol)
    end
    for field_cp in field_cps
        b = sum((VacuumFields.dψ_dZ(coil, field_cp.R, field_cp.Z) * cos(field_cp.θ) + VacuumFields.dψ_dR(coil, field_cp.R, field_cp.Z) * sin(field_cp.θ)) / 2π / field_cp.R for coil in coils)
        @test isapprox(b, field_cp.target; atol=tol)
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
