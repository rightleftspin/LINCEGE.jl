struct ExpansionLattice{C<:Coordinates,N<:Neighbors}
        coordinates::C
        neighbors::N
end

function ExpansionLattice(tiling::Tiling, max_order::Int)
        coordinates = Coordinates(tiling, max_order)
        neighbors = Neighbors(tiling, coordinates)

        ExpansionLattice(
                coordinates,
                neighbors
        )
end

neighbors(expansion_lattice::ExpansionLattice, vs::Vertices) = neighbors(expansion_lattice.neighbors, vs)
all_coordinates(expansion_lattice::ExpansionLattice) = all_coordinates(expansion_lattice.coordinates)
