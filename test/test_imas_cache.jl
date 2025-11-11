using VacuumFields
using Test

@testset "IMAS Coil Caching Infrastructure" begin

    @testset "DerivedCoilData struct" begin
        # Test DerivedCoilData creation
        turns = 100.0
        outline_r = [1.0, 2.0, 2.0, 1.0]
        outline_z = [0.0, 0.0, 1.0, 1.0]

        derived = VacuumFields.DerivedCoilData{Float64}(turns, outline_r, outline_z)

        @test derived.turns_with_sign == 100.0
        @test derived.outline_r == outline_r
        @test derived.outline_z == outline_z
        @test length(derived.outline_r) == 4
        @test length(derived.outline_z) == 4

        # Test with different numeric type
        derived_f32 = VacuumFields.DerivedCoilData{Float32}(Float32(50.0),
                                                             Float32[1.0, 2.0],
                                                             Float32[0.0, 1.0])
        @test derived_f32.turns_with_sign isa Float32
        @test derived_f32.outline_r isa Vector{Float32}
        @test derived_f32.outline_z isa Vector{Float32}
    end

    @testset "CurrentCache struct" begin
        # Test CurrentCache creation
        time = 0.5
        current = 1500.0

        cache = VacuumFields.CurrentCache{Float64}(time, current)

        @test cache.time == 0.5
        @test cache.current_per_turn == 1500.0

        # Test with different input types - should convert to specified type
        cache_from_int = VacuumFields.CurrentCache{Float64}(1, 2000)
        @test cache_from_int.time isa Float64
        @test cache_from_int.current_per_turn isa Float64
        @test cache_from_int.time == 1.0
        @test cache_from_int.current_per_turn == 2000.0

        # Test Float32 cache
        cache_f32 = VacuumFields.CurrentCache{Float32}(Float32(0.5), Float32(1500.0))
        @test cache_f32.time isa Float32
        @test cache_f32.current_per_turn isa Float32
    end

    @testset "DerivedCoilData vector operations" begin
        # Test creating a vector of derived data (as used in actual coils)
        derived_vec = VacuumFields.DerivedCoilData{Float64}[]

        push!(derived_vec, VacuumFields.DerivedCoilData{Float64}(
            100.0, [1.0, 2.0, 2.0, 1.0], [0.0, 0.0, 1.0, 1.0]
        ))
        push!(derived_vec, VacuumFields.DerivedCoilData{Float64}(
            200.0, [1.1, 2.1, 2.1, 1.1], [0.1, 0.1, 1.1, 1.1]
        ))

        @test length(derived_vec) == 2
        @test derived_vec[1].turns_with_sign == 100.0
        @test derived_vec[2].turns_with_sign == 200.0

        # Test accessing cached outline data
        @test length(derived_vec[1].outline_r) == 4
        @test derived_vec[1].outline_r[1] == 1.0
        @test derived_vec[2].outline_z[end] == 1.1
    end

    @testset "CurrentCache time matching semantics" begin
        # Test exact time matching behavior
        time_val = 0.123456789
        current_val = 1234.5678

        cache = VacuumFields.CurrentCache{Float64}(time_val, current_val)

        # Exact match
        @test cache.time == time_val

        # Very close but not exact - should be different
        @test cache.time != (time_val + 1e-15)
        @test cache.time != (time_val - 1e-15)

        # This demonstrates the exact matching semantics needed for get_cached_current
        test_time = 0.5
        cache2 = VacuumFields.CurrentCache{Float64}(test_time, 1000.0)
        @test cache2.time == 0.5  # Exact match
        @test cache2.time != 0.50000001  # Not a match
    end

    @testset "Type conversion in CurrentCache" begin
        # Test that constructor properly converts types

        # Int to Float64
        cache1 = VacuumFields.CurrentCache{Float64}(1, 2000)
        @test cache1.time isa Float64
        @test cache1.current_per_turn isa Float64

        # Float32 to Float64
        cache2 = VacuumFields.CurrentCache{Float64}(Float32(0.5), Float32(1500.0))
        @test cache2.time isa Float64
        @test cache2.current_per_turn isa Float64
        @test cache2.time == 0.5
        @test cache2.current_per_turn == 1500.0

        # Float64 to Float32
        cache3 = VacuumFields.CurrentCache{Float32}(0.5, 1500.0)
        @test cache3.time isa Float32
        @test cache3.current_per_turn isa Float32
    end

    @testset "Cache helper function signatures" begin
        # Verify that cache management functions exist with correct signatures
        @test hasmethod(VacuumFields.has_cached_current, (VacuumFields.GS_IMAS_pf_active__coil,))
        @test hasmethod(VacuumFields.get_cached_current, (VacuumFields.GS_IMAS_pf_active__coil, Real))
        @test hasmethod(VacuumFields.cache_current!, (VacuumFields.GS_IMAS_pf_active__coil, Real, Real))
        @test hasmethod(VacuumFields.clear_current_cache!, (VacuumFields.GS_IMAS_pf_active__coil,))

        @test hasmethod(VacuumFields.has_derived_data, (VacuumFields.GS_IMAS_pf_active__coil,))
        @test hasmethod(VacuumFields.derive_coil_data!, (VacuumFields.GS_IMAS_pf_active__coil,))
        @test hasmethod(VacuumFields.clear_derived_data!, (VacuumFields.GS_IMAS_pf_active__coil,))
    end
end

println("\n✓ All cache infrastructure tests passed!")
println("  - DerivedCoilData struct validated")
println("  - CurrentCache struct validated")
println("  - Type conversion behavior verified")
println("  - Cache management functions exist")
println("\nNote: Integration tests with actual IMAS coils should be added")
println("      to test files that load IMAS data structures.")
