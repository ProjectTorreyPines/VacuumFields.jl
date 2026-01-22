@testset "In-place flux functions" begin

    @testset "ψ! vs ψ" begin
        # Test at multiple points
        test_points = [
            (1.5, 0.0),
            (2.0, 0.5),
            (1.8, -0.3),
            (2.5, 1.0),
        ]

        # Pre-allocate result array
        result = zeros(length(coils))

        for (R, Z) in test_points
            # Compute using in-place version
            VacuumFields.ψ!(result, coils, R, Z)

            # Compute using original allocating version
            expected = [VacuumFields.ψ(coil, R, Z) for coil in coils]

            # Verify identical results
            @test result == expected
            @test all(isfinite, result)
        end
    end

    @testset "dψ_dR! vs dψ_dR" begin
        test_points = [
            (1.5, 0.0),
            (2.0, 0.5),
            (1.8, -0.3),
            (2.5, 1.0),
        ]

        result = zeros(length(coils))

        for (R, Z) in test_points
            # In-place version
            VacuumFields.dψ_dR!(result, coils, R, Z)

            # Original version
            expected = [VacuumFields.dψ_dR(coil, R, Z) for coil in coils]

            @test result == expected
            @test all(isfinite, result)
        end
    end

    @testset "dψ_dZ! vs dψ_dZ" begin
        test_points = [
            (1.5, 0.0),
            (2.0, 0.5),
            (1.8, -0.3),
            (2.5, 1.0),
        ]

        result = zeros(length(coils))

        for (R, Z) in test_points
            # In-place version
            VacuumFields.dψ_dZ!(result, coils, R, Z)

            # Original version
            expected = [VacuumFields.dψ_dZ(coil, R, Z) for coil in coils]

            @test result == expected
            @test all(isfinite, result)
        end
    end

    @testset "d2ψ_dZ2! vs d2ψ_dZ2" begin
        test_points = [
            (1.5, 0.0),
            (2.0, 0.5),
            (1.8, -0.3),
            (2.5, 1.0),
        ]

        result = zeros(length(coils))

        for (R, Z) in test_points
            # In-place version
            VacuumFields.d2ψ_dZ2!(result, coils, R, Z)

            # Original version
            expected = [VacuumFields.d2ψ_dZ2(coil, R, Z) for coil in coils]

            @test result == expected
            @test all(isfinite, result)
        end
    end

    @testset "flux_on_grid! vs flux_on_grid" begin
        # Create a test grid
        Rs = range(0.8, 2.8; length=25)
        Zs = range(-1.5, 1.5; length=25)

        # Compute Green's function table
        Gtable = VacuumFields.Green_table(Rs, Zs, coils)

        # Allocating version
        Ψ_alloc = VacuumFields.flux_on_grid(Gtable, Rs, Zs, coils)

        # In-place version
        Ψ_inplace = zeros(eltype(Ψ_alloc), size(Ψ_alloc))
        VacuumFields.flux_on_grid!(Ψ_inplace, Gtable, Rs, Zs, coils)

        # Verify identical results
        @test Ψ_inplace == Ψ_alloc
        @test size(Ψ_inplace) == (length(Rs), length(Zs))
        @test all(isfinite, Ψ_inplace)
    end

    @testset "MultiCoil support" begin
        # Create a MultiCoil
        orientation = ones(Int, length(coils))
        orientation[1:2:end] .= -1
        mcoil = VacuumFields.MultiCoil(coils, orientation)
        VacuumFields.set_current_per_turn!(mcoil, Icpt)

        R, Z = 2.0, 0.0

        # Test ψ! with vector of MultiCoil
        mcoils = [mcoil]
        result = zeros(length(mcoils))
        VacuumFields.ψ!(result, mcoils, R, Z)

        # Compare to original
        expected = VacuumFields.ψ(mcoil, R, Z)
        @test result[1] == expected
    end
end

@testset "In-place allocation tests" begin
    # Optional: verify that in-place versions don't allocate
    # This is useful for performance validation

    R, Z = 2.0, 0.0
    result = zeros(length(coils))

    # Warm-up run to compile
    VacuumFields.ψ!(result, coils, R, Z)

    # Check allocations (should be minimal or zero)
    allocs = @allocated VacuumFields.ψ!(result, coils, R, Z)
    @test allocs == 0 || allocs < 100  # Allow minimal allocations

    # Compare to allocating version
    allocs_original = @allocated [VacuumFields.ψ(coil, R, Z) for coil in coils]
    @test allocs < allocs_original  # In-place should allocate less
end

@testset "In-place consistency with summed values" begin
    # Test that summing in-place results matches direct sum computation
    R, Z = 2.0, 0.0

    # original allocating version
    result_ori = [VacuumFields.ψ(coil, R, Z) for coil in coils]
    sum_ori = sum(result_ori)

    # Using direct sum (as in MultiCoil)
    sum_direct = sum(VacuumFields.ψ(coil, R, Z) for coil in coils) 

    # Using in-place functions
    result = zeros(length(coils))
    VacuumFields.ψ!(result, coils, R, Z)
    sum_inplace = sum(result)


    @test sum_ori == sum_inplace # exactly equal
    
    # Direct sum may differ slightly (truncation error) due to different internal order of operations
    @test isapprox(sum_ori, sum_direct, rtol=1e-12) 
end

println("\n✓ All in-place function tests passed!")
println("  In-place versions produce identical results to original implementations")
