using VacuumFields
using Test
import VacuumFields: IMAS
import VacuumFields.IMAS.ForwardDiff: Dual, value, partials

# Coil Green's functions and fluxes with dual numbers, against central finite differences

@testset "coils with dual numbers" begin
    seed(v, k) = Dual{Nothing}(v, ntuple(i -> Float64(i == k), 2)...)
    R0, Z0, w, hgt = 1.5, 0.4, 0.2, 0.3
    function corners(Rc, Zc) # both coordinates in the same number type
        T = promote_type(typeof(Rc), typeof(Zc))
        return (T[Rc - w / 2, Rc + w / 2, Rc + w / 2, Rc - w / 2], T[Zc - hgt / 2, Zc - hgt / 2, Zc + hgt / 2, Zc + hgt / 2])
    end

    # derivative with respect to (1) the radius of the coil and (2) the radius of the observation point
    G(Rc, R) = VacuumFields.Green(VacuumFields.QuadCoil(corners(Rc, Z0)..., 1.0), R, -0.2)
    g = G(seed(R0, 1), seed(2.3, 2))
    @test value(g) ≈ G(R0, 2.3)
    h = 1e-5
    @test partials(g, 1) ≈ (G(R0 + h, 2.3) - G(R0 - h, 2.3)) / (2h) rtol = 1e-6
    @test partials(g, 2) ≈ (G(R0, 2.3 + h) - G(R0, 2.3 - h)) / (2h) rtol = 1e-6

    # a coil whose geometry carries derivatives and whose current does not (and the reverse)
    geom = VacuumFields.QuadCoil(seed.(corners(R0, Z0)[1], 1), seed.(corners(R0, Z0)[2], 2), 1e3)
    @test VacuumFields.Green(geom, 2.3, -0.2) isa Dual
    curr = VacuumFields.QuadCoil(corners(R0, Z0)..., seed(1e3, 1))
    @test VacuumFields.Green(curr, 2.3, -0.2) ≈ G(R0, 2.3)

    # coils in series: the common current takes the number type of the coils
    both = VacuumFields.QuadCoil(seed.(corners(R0, Z0)[1], 1), seed.(corners(R0, Z0)[2], 2), seed(1e3, 1))
    multi = VacuumFields.MultiCoil([both, deepcopy(both)], [1, -1])
    @test VacuumFields.current_per_turn(multi.coils[1]) isa Real
    VacuumFields.set_current_per_turn!(multi, seed(2e3, 1))
    @test value(VacuumFields.current_per_turn(multi.coils[1])) ≈ 2e3
end
