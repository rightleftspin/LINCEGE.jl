struct Lattice{E<:ExpansionLattice,R<:RealSpaceLattice,M<:AbstractAdjacencyMatrix}
        expansion_lattice::E
        real_space_lattice::R
        masking_matrix::M
end

function Lattice(tiling::Tiling, max_order::Int)
        expansion_lattice = ExpansionLattice(tiling, max_order)
        real_space_lattice = RealSpaceLattice(tiling, expansion_lattice)
        masking_matrix = MaskingAdjacencyMatrix(expansion_lattice, real_space_lattice)

        Lattice(
                expansion_lattice,
                real_space_lattice,
                masking_matrix,
        )
end
