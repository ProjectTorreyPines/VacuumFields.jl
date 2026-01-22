using VacuumFields
using Test

@testset "IMAS Coil Caching Infrastructure" begin

    @testset "ElementCache struct" begin
        # Test ElementCache creation
        turns = 100.0
        outline_r = [1.0, 2.0, 2.0, 1.0]
        outline_z = [0.0, 0.0, 1.0, 1.0]
        geom = (r=1.5, z=0.5, width=1.0, height=1.0)

        elem_cache = VacuumFields.ElementCache{Float64}(turns, outline_r, outline_z, geom)

        @test elem_cache.turns_with_sign == 100.0
        @test elem_cache.outline_r == outline_r
        @test elem_cache.outline_z == outline_z
        @test elem_cache.source_geometry == geom
        @test length(elem_cache.outline_r) == 4
        @test length(elem_cache.outline_z) == 4

        # Test with different numeric type
        geom_f32 = (r=Float32(1.5), z=Float32(0.5), width=Float32(1.0), height=Float32(1.0))
        elem_cache_f32 = VacuumFields.ElementCache{Float32}(Float32(50.0),
                                                             Float32[1.0, 2.0],
                                                             Float32[0.0, 1.0],
                                                             geom_f32)
        @test elem_cache_f32.turns_with_sign isa Float32
        @test elem_cache_f32.outline_r isa Vector{Float32}
        @test elem_cache_f32.outline_z isa Vector{Float32}
        @test elem_cache_f32.source_geometry isa NamedTuple
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

    @testset "ElementCache vector operations" begin
        # Test creating a vector of element caches (as used in actual coils)
        elements_cache = VacuumFields.ElementCache{Float64}[]

        push!(elements_cache, VacuumFields.ElementCache{Float64}(
            100.0, [1.0, 2.0, 2.0, 1.0], [0.0, 0.0, 1.0, 1.0],
            (r=1.5, z=0.5, width=1.0, height=1.0)
        ))
        push!(elements_cache, VacuumFields.ElementCache{Float64}(
            200.0, [1.1, 2.1, 2.1, 1.1], [0.1, 0.1, 1.1, 1.1],
            (r=1.6, z=0.6, width=1.0, height=1.0)
        ))

        @test length(elements_cache) == 2
        @test elements_cache[1].turns_with_sign == 100.0
        @test elements_cache[2].turns_with_sign == 200.0

        # Test accessing cached outline data
        @test length(elements_cache[1].outline_r) == 4
        @test elements_cache[1].outline_r[1] == 1.0
        @test elements_cache[2].outline_z[end] == 1.1

        # Test source geometry
        @test elements_cache[1].source_geometry.r == 1.5
        @test elements_cache[2].source_geometry.width == 1.0
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
        # is_cache_valid has three methods via multiple dispatch
        @test hasmethod(VacuumFields.is_cache_valid, (Vector{VacuumFields.ElementCache{Float64}}, Any))
        @test hasmethod(VacuumFields.is_cache_valid, (VacuumFields.CurrentCache, Real))

        # Cache management functions
        @test hasmethod(VacuumFields.cache_current!, (VacuumFields.GS_IMAS_pf_active__coil, Real, Real))
        @test hasmethod(VacuumFields.ensure_valid_elements_cache!, (VacuumFields.GS_IMAS_pf_active__coil,))
        @test hasmethod(VacuumFields.clear_elements_cache!, (VacuumFields.GS_IMAS_pf_active__coil,))
    end

    @testset "is_cache_valid behavior" begin
        # Test is_cache_valid with CurrentCache
        cache_t1 = VacuumFields.CurrentCache{Float64}(0.5, 1500.0)
        cache_t2 = VacuumFields.CurrentCache{Float64}(1.0, 2000.0)
        cache_nan = VacuumFields.CurrentCache{Float64}(NaN, NaN)

        # Test exact time matching semantics
        @test VacuumFields.is_cache_valid(cache_t1, 0.5)
        @test !VacuumFields.is_cache_valid(cache_t1, 0.5000001)
        @test !VacuumFields.is_cache_valid(cache_t2, 0.5)

        # Test NaN cache (invalid)
        @test !VacuumFields.is_cache_valid(cache_nan, 0.5)
        @test !VacuumFields.is_cache_valid(cache_nan, NaN)
    end
end

println("\n✓ All cache infrastructure tests passed!")
println("  - ElementCache struct validated")
println("  - CurrentCache struct validated")
println("  - Type conversion behavior verified")
println("  - Cache management functions exist")
println("  - is_cache_valid direct cache validation tested")
println("\nNote: Integration tests with actual IMAS coils should be added")
println("      to test files that load IMAS data structures.")
