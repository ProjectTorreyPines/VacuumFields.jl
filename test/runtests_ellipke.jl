using VacuumFields
using Test

@testset "Elliptic Integrals Correctness" begin

    # Test values across different ranges and edge cases
    test_values_positive = [
        0.0,    # Edge: zero
        0.05,   # Range 1: 0.0-0.1
        0.15,   # Range 2: 0.1-0.2
        0.25,   # Range 3: 0.2-0.3
        0.35,   # Range 4: 0.3-0.4
        0.45,   # Range 5: 0.4-0.5
        0.55,   # Range 6: 0.5-0.6
        0.65,   # Range 7: 0.6-0.7
        0.75,   # Range 8: 0.7-0.8
        0.825,  # Range 9: 0.8-0.85
        0.875,  # Range 10: 0.85-0.9
        0.95,   # Special: >= 0.9
        0.99,   # Near singularity
        0.999,  # Very near singularity
    ]

    # Negative values chosen so that x = m/(m-1) maps to each polynomial table
    test_values_negative = [
        -0.05,   # x ≈ 0.048  → Table 2 (0.0-0.1)
        -0.17,   # x ≈ 0.145  → Table 3 (0.1-0.2)
        -0.33,   # x ≈ 0.248  → Table 4 (0.2-0.3)
        -0.54,   # x ≈ 0.351  → Table 5 (0.3-0.4)
        -0.82,   # x ≈ 0.451  → Table 6 (0.4-0.5)
        -1.22,   # x ≈ 0.550  → Table 7 (0.5-0.6)
        -1.86,   # x ≈ 0.650  → Table 8 (0.6-0.7)
        -3.0,    # x = 0.75   → Table 9 (0.7-0.8)
        -4.7,    # x ≈ 0.825  → Table 10 (0.8-0.85)
        -7.0,    # x = 0.875  → Table 11 (0.85-0.9)
        -19.0,   # x = 0.95   → Special (0.9-1.0)
    ]

    test_values_all = vcat(test_values_negative, test_values_positive)

    @testset "ellipke vs ellipk/ellipe (all m)" begin
        for m in test_values_all
            # Original individual functions
            Km_orig = VacuumFields.ellipk(m)
            Em_orig = VacuumFields.ellipe(m)

            # Combined function
            Km_new, Em_new = VacuumFields.ellipke(m)

            # Test exact bit-level equality
            @test Km_new == Km_orig
            @test Em_new == Em_orig
        end
    end

    @testset "ellipke_nonneg correctness (positive m only)" begin
        for m in test_values_positive
            # Test against ellipk/ellipe
            Km_ref = VacuumFields.ellipk(m)
            Em_ref = VacuumFields.ellipe(m)

            Km_test, Em_test = VacuumFields.ellipke_nonneg(m)

            @test Km_test == Km_ref
            @test Em_test == Em_ref
        end
    end

    @testset "ellipke correctness (negative m only)" begin
        for m in test_values_negative
            # Test against ellipk/ellipe - bit-level equality
            Km_ref = VacuumFields.ellipk(m)
            Em_ref = VacuumFields.ellipe(m)

            Km_test, Em_test = VacuumFields.ellipke(m)

            @test Km_test == Km_ref
            @test Em_test == Em_ref
        end
    end

    @testset "ellipke negative m handling" begin
        for m in test_values_negative
            # Should handle negative m correctly
            Km, Em = VacuumFields.ellipke(m)

            # Basic sanity checks
            @test isfinite(Km)
            @test isfinite(Em)
            @test Km > 0  # K(m) is always positive
            @test Em > 0  # E(m) is always positive for negative m
        end
    end

    @testset "Negative m transformation to polynomial tables" begin
        # For negative m, the transformation is: x = m / (m - 1)
        # Test values chosen so that transformed x hits each polynomial table

        test_cases = [
            # (m, expected_table_range, description)
            (-0.05,   "Table 2 (0.0-0.1)",    "x ≈ 0.048"),
            (-0.17,   "Table 3 (0.1-0.2)",    "x ≈ 0.145"),
            (-0.33,   "Table 4 (0.2-0.3)",    "x ≈ 0.248"),
            (-0.54,   "Table 5 (0.3-0.4)",    "x ≈ 0.351"),
            (-0.82,   "Table 6 (0.4-0.5)",    "x ≈ 0.451"),
            (-1.22,   "Table 7 (0.5-0.6)",    "x ≈ 0.550"),
            (-1.86,   "Table 8 (0.6-0.7)",    "x ≈ 0.650"),
            (-3.0,    "Table 9 (0.7-0.8)",    "x = 0.75"),
            (-4.7,    "Table 10 (0.8-0.85)",  "x ≈ 0.825"),
            (-7.0,    "Table 11 (0.85-0.9)",  "x = 0.875"),
            (-19.0,   "Special (0.9-1.0)",    "x = 0.95"),
        ]

        for (m, table_desc, x_desc) in test_cases
            # Verify transformation
            x_transformed = m / (m - 1)

            # Test against reference implementation
            Km_ref = VacuumFields.ellipk(m)
            Em_ref = VacuumFields.ellipe(m)

            Km_test, Em_test = VacuumFields.ellipke(m)

            @test Km_test == Km_ref
            @test Em_test == Em_ref
            @test isfinite(Km_test)
            @test isfinite(Em_test)
            @test Km_test > 0
            @test Em_test > 0
        end
    end

    @testset "ellipke edge cases" begin
        # Test m = 0
        Km, Em = VacuumFields.ellipke(0.0)
        @test Km ≈ π/2
        @test Em ≈ π/2

        # Test m = 1 (should return Inf, 1)
        Km, Em = VacuumFields.ellipke(1.0)
        @test isinf(Km)
        @test Em ≈ 1.0

        # Test m very close to 1
        Km, Em = VacuumFields.ellipke(0.9999)
        @test isfinite(Km)
        @test Km > π/2  # Should be larger than K(0)
        @test isfinite(Em)

        # Test m = 1.0 + eps() (should still work, treated as ≈ 1)
        Km, Em = VacuumFields.ellipke(1.0 + eps())
        @test isinf(Km)
        @test Em ≈ 1.0

        # Test m = nextfloat(1.0) (should still work, treated as ≈ 1)
        Km, Em = VacuumFields.ellipke(nextfloat(1.0))
        @test isinf(Km)
        @test Em ≈ 1.0
    end

    @testset "Domain errors for m > 1" begin
        # Values clearly beyond the valid domain should throw DomainError
        @test_throws DomainError VacuumFields.ellipke_nonneg(1.01)
        @test_throws DomainError VacuumFields.ellipke_nonneg(1.1)
        @test_throws DomainError VacuumFields.ellipke_nonneg(2.0)
        @test_throws DomainError VacuumFields.ellipke_nonneg(10.0)

        # But values very close to 1 (within floating point tolerance) should work
        # because x ≈ 1.0 is checked
        Km, Em = VacuumFields.ellipke_nonneg(1.0 + eps())
        @test isinf(Km)
        @test Em ≈ 1.0
    end

    @testset "Consistency across ranges" begin
        # Test boundary values between polynomial ranges
        boundaries = [
            0.099, 0.1, 0.101,    # Table 2-3 boundary
            0.199, 0.2, 0.201,    # Table 3-4 boundary
            0.299, 0.3, 0.301,    # Table 4-5 boundary
            0.399, 0.4, 0.401,    # Table 5-6 boundary
            0.499, 0.5, 0.501,    # Table 6-7 boundary
            0.599, 0.6, 0.601,    # Table 7-8 boundary
            0.699, 0.7, 0.701,    # Table 8-9 boundary
            0.799, 0.8, 0.801,    # Table 9-10 boundary
            0.849, 0.85, 0.851,   # Table 10-11 boundary
            0.899, 0.9, 0.901,    # Table 11-special boundary
        ]

        for m in boundaries
            Km_ref = VacuumFields.ellipk(m)
            Em_ref = VacuumFields.ellipe(m)

            Km_test, Em_test = VacuumFields.ellipke(m)

            @test Km_test == Km_ref
            @test Em_test == Em_ref
        end
    end
end
