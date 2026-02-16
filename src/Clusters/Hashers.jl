struct TranslationHasher <: AbstractHasher
    hashing_matrix::AbstractMatrix{Int}
end

function TranslationHasher(lattice::SiteExpansionLattice)
    pwd_matrix = pairwise_direction(get_coordinates(lattice))
    hashing_matrix::Matrix{Int}, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
    diag_matrix = diagm(get_labels(lattice) .+ max_dir)
    TranslationHasher(2 .^ (hashing_matrix + diag_matrix))
end

ghash(h::TranslationHasher, evs::ExpansionVertices) = hash(sum(h.hashing_matrix[evs, evs], dims=2))

struct IsomorphicHasher <: AbstractHasher
    hashing_matrix::AbstractMatrix{Int}
    labels::AbstractVector{Int}
    is_weighted::Bool
end

function IsomorphicHasher(lattice::SiteExpansionLattice)
    hashing_matrix = bond_matrix(lattice)
    if length(unique(hashing_matrix)) > 2
        IsomorphicHasher(hashing_matrix, get_site_colors(lattice), true)
    end

    IsomorphicHasher(hashing_matrix, get_site_colors(lattice), false)
end

function ghash(hasher::IsomorphicHasher, evs::ExpansionVertices)
    if hasher.is_weighted
        (h, p) = weighted_iso_hash(
            hasher.hashing_matrix[evs, evs],
            hasher.labels[evs]
        )
    else
        (h, p) = unweighted_iso_hash(
            hasher.hashing_matrix[evs, evs],
            hasher.labels[evs]
        )
    end
    h
end

"""
Finds the canonical ordering of an edge labeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function weighted_iso_hash(weighted_adj_mat::AbstractMatrix{Int}, labels::AbstractVector{Int})

    # This is to get edge-labels to work
    nauty_labels = vcat(labels, zeros(Int64, fld(count(>(1), weighted_adj_mat), 2)))

    unweighted_adjacency_matrix = zeros(Int64, length(nauty_labels), length(nauty_labels))

    current_aux_vert = size(weighted_adj_mat, 1) + 1
    for j = 1:size(weighted_adj_mat, 2)
        for i = j:size(weighted_adj_mat, 1)

            if (@inbounds(weighted_adj_mat[i, j]) == 1)
                # Add an edge if there is already an edge
                @inbounds unweighted_adjacency_matrix[i, j] = 1
                @inbounds unweighted_adjacency_matrix[j, i] = 1
            elseif (@inbounds(weighted_adj_mat[i, j]) != 0)
                # Add an edge to the auxilary vertex here
                @inbounds unweighted_adjacency_matrix[i, current_aux_vert] = 1
                @inbounds unweighted_adjacency_matrix[j, current_aux_vert] = 1
                @inbounds unweighted_adjacency_matrix[current_aux_vert, i] = 1
                @inbounds unweighted_adjacency_matrix[current_aux_vert, j] = 1
                # Color the aux vertex the same as the edge
                @inbounds nauty_labels[current_aux_vert] = weighted_adj_mat[i, j]
                current_aux_vert += 1
            end
        end
    end

    nauty_graph = NautyGraph(unweighted_adjacency_matrix, vertex_labels=nauty_labels)
    # Canonize and find the corresponding permutation
    permutation = canonize!(nauty_graph)

    # TODO: Add cluster symmetries correctly

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes from the original graph to the
    # canonized graph
    return (NautyGraphs.ghash(nauty_graph), @inbounds(permutation[1:size(weighted_adj_mat, 1)]))

end

"""
Finds the canonical ordering of an edge unlabeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function unweighted_iso_hash(unweighted_adj_mat::AbstractMatrix{Int}, labels::AbstractVector{Int})

    nauty_graph = NautyGraph(unweighted_adj_mat, vertex_labels=labels)
    # Canonize and find the corresponding permutation
    permutation = canonize!(nauty_graph)

    # TODO: Add cluster symmetries correctly

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes from the original graph to the
    # canonized graph
    return (NautyGraphs.ghash(nauty_graph), @inbounds(permutation[1:size(unweighted_adj_mat, 1)]))

end
