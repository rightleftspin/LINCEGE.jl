begin # Hashing Functions

    """
    Takes a cluster and finds the translationally invariant hash of it. This assumes that each
    site in the fundamental tiling unit is colored differently. In this way, the coloring of
    the translational units supersedes the labeling of each site.

    Please note that this function requires all the vertices to be sorted in dimensional order,
    ie. x, y, z, ...
    """
    function translational_pruning(cluster::Cluster)
        (
            hash(
                vec(
                    sum(weight -> 2^weight, direction_adjacency_matrix(cluster), dims = 2),
                ) + (1 .// (all_translation_labels(cluster) .+ 1)),
            ),
            nothing,
        )
    end

    """
    Finds the canonical ordering of an edge labeled cluster using Nauty, this function
    returns the permutation to rearrange the cluster used by Nauty
    """
    function isomorphic_pruning(cluster::Cluster{<:Any,true})

        weighted_adj_mat = weighted_adjacency_matrix(cluster)
        # This is to get edge-labels to work
        nauty_labels = vcat(
            all_vertex_labels(cluster),
            zeros(Int64, fld(count(>(1), weighted_adj_mat), 2)),
        )

        unweighted_adjacency_matrix =
            zeros(Int64, length(nauty_labels), length(nauty_labels))

        current_aux_vert = nv(cluster) + 1
        for j = 1:size(weighted_adj_mat, 2)
            for i = j:size(weighted_adj_mat, 1)

                if (weighted_adj_mat[i, j] == 1)
                    # Add an edge if there is already an edge
                    unweighted_adjacency_matrix[i, j] = 1
                    unweighted_adjacency_matrix[j, i] = 1
                elseif (weighted_adj_mat[i, j] != 0)
                    # Add an edge to the auxilary vertex here
                    unweighted_adjacency_matrix[i, current_aux_vert] = 1
                    unweighted_adjacency_matrix[j, current_aux_vert] = 1
                    unweighted_adjacency_matrix[current_aux_vert, i] = 1
                    unweighted_adjacency_matrix[current_aux_vert, j] = 1
                    # Color the aux vertex the same as the edge
                    nauty_labels[current_aux_vert] = weighted_adj_mat[i, j]
                    current_aux_vert += 1
                end
            end
        end

        nauty_graph = NautyGraph(unweighted_adjacency_matrix, nauty_labels)
        # Canonize and find the corresponding permutation
        #permutation = canonize!(nauty_graph)
        permutation = canonical_permutation(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), permutation[1:nv(cluster)])

    end

    """
    Finds the canonical ordering of an edge unlabeled cluster using Nauty, this function
    returns the permutation to rearrange the cluster used by Nauty
    """
    function isomorphic_pruning(cluster::Cluster{<:Any,false})

        nauty_graph =
            NautyGraph(weighted_adjacency_matrix(cluster), all_vertex_labels(cluster))
        # Canonize and find the corresponding permutation
        permutation = canonize!(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), permutation[1:nv(cluster)])

    end

    # TODO: Need to write the symmetric pruning function

end
