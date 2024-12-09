
using NLCE, Test

using Base.Iterators

@testset verbose = true "NLCE" begin

    # A few different lattices to test on
    square_lattice = Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[1, 0], [0, 1]])

    triangular_lattice =
        Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[0.5, sqrt(3) / 2], [1, 0]])
    
    kagome_lattice =
        Dict("Basis" => [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]], "Primitive Vectors" => [[2, 0], [-1, sqrt(3)]])

    @testset "NLCELattice" begin
        # Here I will test a variety of lattices at nearest neighbor up till order 6
        # to see if the NLCELattice struct works. Most issues usually present themselves
        # by order 5.

        neighborhood = [1]
        max_order = 6

        square_nlce_lattice = NLCE.NLCELattice(
            square_lattice["Basis"],
            square_lattice["Primitive Vectors"],
            neighborhood,
            max_order,
        )

        triangular_nlce_lattice = NLCE.NLCELattice(
            triangular_lattice["Basis"],
            triangular_lattice["Primitive Vectors"],
            neighborhood,
            max_order,
        )

        kagome_nlce_lattice = NLCE.NLCELattice(
            kagome_lattice["Basis"],
            kagome_lattice["Primitive Vectors"],
            neighborhood,
            max_order,
        )

        number_verts = (2 * max_order + 1)^2

        @test NLCE.nv(square_nlce_lattice) == number_verts
        @test NLCE.nv(triangular_nlce_lattice) == number_verts

        @test NLCE.center(square_nlce_lattice) == [fld(number_verts, 2)]
        @test NLCE.center(triangular_nlce_lattice) == [fld(number_verts, 2)]

        @test length(NLCE.neighbors(square_nlce_lattice, fld(number_verts, 2))) == 4
        @test length(NLCE.neighbors(triangular_nlce_lattice, fld(number_verts, 2))) == 6

        @test length(NLCE.neighbors(square_nlce_lattice, [fld(number_verts, 2), 1])) == 2
        @test length(NLCE.neighbors(triangular_nlce_lattice, [fld(number_verts, 2), 1])) ==
              2

        @test sum(length, NLCE.neighbors(square_nlce_lattice, [fld(number_verts, 2), 1])) ==
              6
        @test sum(
            length,
            NLCE.neighbors(triangular_nlce_lattice, [fld(number_verts, 2), 1]),
        ) == 8

    end
end
