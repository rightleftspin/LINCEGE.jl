struct RealSpaceLattice{C<:Coordinates,D<:AbstractAdjacencyMatrix,B<:AbstractAdjacencyMatrix}
        coordinates::C
        direction_matrix::D
        bond_matrix::B
end

function RealSpaceLattice(tiling::Tiling, expansion_lattice::ExpansionLattice)
        coordinates = Coordinates(tiling, all_coordinates(expansion_lattice))
        direction_matrix = DirectionMatrix(tiling, coordinates)
        bond_matrix = BondMatrix(tiling, coordinates)

        RealSpaceLattice(
                coordinates,
                direction_matrix,
                bond_matrix,
        )
end
