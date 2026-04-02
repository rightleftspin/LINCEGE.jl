struct WeakClusterExpansionLattice <: AbstractClusterExpansionLattice
        max_order::UInt8

        expansion_unit_cell::ExpansionUnitCell
        centers::ExpansionVertices
        neighbor_list::Vector{ExpansionVertices{Int}}

        lattice_coordinates::Matrix{Int}
        translation_labels::Vector{Int}
        site_colors::Vector{Int}
        adj_matrix::Matrix{Int}

        connections::WeakClusterConnections
end

function WeakClusterExpansionLattice(max_order::Int, expansion_unit_cell::ExpansionUnitCell)
        @assert max_order > 0 "max_order must be a positive integer"

        expansion_coordinates = generate_coordinates(max_order, length(basis_size(expansion_unit_cell)), dimension(expansion_unit_cell))
        neighbor_list = generate_neighbor_list(expansion_coordinates, expansion_unit_cell)
        ctrs = ExpansionVertices(find_centers(expansion_coordinates))

        lattice_coordinates = []
        for i in 1:length(basis_size(expansion_unit_cell))
                coords = generate_coordinates(max_order, basis_size(expansion_unit_cell)[i], dimension(expansion_unit_cell))
                push!(lattice_coordinates, vcat(coords[1:dimension(expansion_unit_cell), :], ones(Int, size(coords, 2))' * i, coords[end, :]'))
        end
        lattice_coordinates = hcat(lattice_coordinates...)

        connections_vec, rev_connections, masking_matrix, unique_inds = generate_weak_connections(expansion_coordinates, lattice_coordinates, expansion_unit_cell)

        adj_matrix = generate_adj_matrix_weak(lattice_coordinates, expansion_unit_cell, unique_inds)

        lattice_coordinates = lattice_coordinates[:, unique_inds]
        translation_labels = [expansion_unit_cell.translation_labels[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        site_colors = [expansion_unit_cell.site_colors[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]

        return WeakClusterExpansionLattice(
                UInt8(max_order),
                expansion_unit_cell,
                ctrs,
                neighbor_list,
                lattice_coordinates,
                translation_labels,
                site_colors,
                adj_matrix,
                WeakClusterConnections(connections_vec, rev_connections, masking_matrix, basis_size(expansion_unit_cell))
        )
end

centers(lattice::WeakClusterExpansionLattice) = lattice.centers
max_order(lattice::WeakClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::WeakClusterExpansionLattice) = length(unique(Iterators.flatten(lattice.expansion_unit_cell.translation_labels)))
n_site_colors(lattice::WeakClusterExpansionLattice) = length(unique(Iterators.flatten(lattice.expansion_unit_cell.site_colors)))
neighbors(lattice::WeakClusterExpansionLattice, vs::ExpansionVertices) = union(ExpansionVertices(), lattice.neighbor_list[vs])
get_coordinates(lattice::WeakClusterExpansionLattice) = shift_unit_cell(lattice.expansion_unit_cell, lattice.lattice_coordinates)
get_labels(lattice::WeakClusterExpansionLattice) = lattice.translation_labels
get_site_colors(lattice::WeakClusterExpansionLattice) = lattice.site_colors
connections(lattice::WeakClusterExpansionLattice) = lattice.connections
bond_matrix(lattice::WeakClusterExpansionLattice) = lattice.adj_matrix
