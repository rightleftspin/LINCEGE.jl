struct StrongClusterExpansionLattice <: AbstractClusterExpansionLattice
        max_order::UInt8
        n_unique_sites::UInt8
        n_site_colors::UInt8
        expansion_unit_cell::UnitCell
        expansion_coordinates::Matrix{Int}
        lattice_unit_cells::Vector{UnitCell}
        lattice_coordinates::Vector{Matrix{Int}}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{ExpansionVertices{Int}}
        connections::StrongClusterConnections
end

function StrongClusterExpansionLattice(max_order::Int, expansion_unit_cell::UnitCell, lattice_unit_cells::Vector{UnitCell})
        @assert max_order > 0 "max_order must be a positive integer"
        @assert length(lattice_unit_cells) == basis_size(expansion_unit_cell) "Need one lattice_unit_cell per expansion site type"

        return StrongClusterExpansionLattice()
end

centers(lattice::StrongClusterExpansionLattice) = ExpansionVertices(find_centers(lattice.expansion_coordinates))
max_order(lattice::StrongClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::StrongClusterExpansionLattice) = lattice.n_unique_sites
n_site_colors(lattice::StrongClusterExpansionLattice) = lattice.n_site_colors
neighbors(lattice::StrongClusterExpansionLattice, vs::ExpansionVertices) = union(ExpansionVertices(), lattice.neighbor_list[vs])
get_coordinates(lattice::StrongClusterExpansionLattice) = hcat([shift_unit_cell(uc, cs) for (uc, cs) in zip(lattice.lattice_unit_cells, lattice.lattice_coordinates)]...)
get_labels(lattice::StrongClusterExpansionLattice) = vcat([cs[end, :] for cs in lattice.lattice_coordinates]...)
get_site_colors(lattice::StrongClusterExpansionLattice) = vcat([uc.site_colors[cs[end, :]] for (uc, cs) in zip(lattice.lattice_unit_cells, lattice.lattice_coordinates)]...)
bond_matrix(lattice::StrongClusterExpansionLattice) = lattice.adj_matrix
connections(lattice::StrongClusterExpansionLattice) = lattice.connections
