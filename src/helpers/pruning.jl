# Need to do
#
#       isomorphic_pruning:
#                       - account for edge weights

"""
These are pruning functions to be used inside step 2 of the algorithm,
these functions will always take in a cluster in some form and return
a hash of the cluster and a potential permutation of the cluster. If
the function returns nothing, this means that no permutation is necessary.
"""

"""
Takes a cluster and finds the translationally invariant
form of it. This does not return a hash.

Inputs: 
      cluster: takes a cluster struct, uses the direction
      weights matrix in the cluster

Output:
      Vector of cluster translational form
"""
function translational_form(cluster::AbstractNLCECluster)
    cluster_perm = sortperm(vertices(cluster))
    (vec(
        sum(
            weight -> 2^weight,
            direction_matrix(cluster)[cluster_perm, :],
            dims = 2,
        ),
    ),
     cluster_perm)
end

"""
Takes a cluster and finds the translationally invariant hash 
of it.

Inputs: 
      cluster: takes a cluster struct, uses the direction
      weights matrix in the cluster

Output:
      Tuple of cluster hash and nothing
"""
function translational_pruning(cluster::AbstractNLCECluster)
    form, perm = translational_form(cluster)
    (hash(form + (1 .// (labels(cluster)[perm] .+ 1))), nothing)

end

"""
Finds the canonical ordering of the cluster using Nauty, this
function rearranges the cluster according to the permutation
used by Nauty

Inputs: 
      cluster: takes in a cluster struct, uses the
      edge_weighted_matrix inside it

Output:
      Tuple of cluster hash and permutation from nauty
"""
function isomorphic_pruning(cluster::AbstractNLCECluster)

    # Create an empty DenseNautyGraph

    # This is all to get weighted edges to work
    nauty_labels = vcat(label(cluster, vertices(cluster)),
                        zeros(Int64, fld(count(>(1), edge_weighted_matrix(cluster)), 2)))

    unweighted_adjacency_matrix = zeros(Int64, length(nauty_labels), length(nauty_labels))

    current_aux_vert = nv(cluster) + 1
    for j = 1:size(edge_weighted_matrix(cluster),2)
        for i = j:size(edge_weighted_matrix(cluster),1)

            if (edge_weighted_matrix(cluster)[i, j] == 1)
            # Add an edge if there is already an edge
                unweighted_adjacency_matrix[i, j] = 1
                unweighted_adjacency_matrix[j, i] = 1
            elseif (edge_weighted_matrix(cluster)[i, j] !=0)
            # Add an edge to the auxilary vertex here
                unweighted_adjacency_matrix[i, current_aux_vert] = 1
                unweighted_adjacency_matrix[j, current_aux_vert] = 1
                unweighted_adjacency_matrix[current_aux_vert, i] = 1
                unweighted_adjacency_matrix[current_aux_vert, j] = 1
            # Color the aux vertex the same as the edge
                nauty_labels[current_aux_vert] = edge_weighted_matrix(cluster)[i, j]
                current_aux_vert += 1
            end
        end
    end

    #nauty_graph =
    #   NautyGraph(edge_weighted_matrix(cluster), label(cluster, vertices(cluster)))

   # println(unweighted_adjacency_matrix)
   # println(nauty_labels)
    nauty_graph =
       NautyGraph(unweighted_adjacency_matrix, nauty_labels)
    # Canonize and find the corresponding permutation 
    permutation, _ = canonize!(nauty_graph)
    
    #add_orbits!(cluster, NautyGraphs.orbits(nauty_graph)) 
    add_orbits!(cluster, vertices(cluster)) 

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes from the original graph to the
    # canonized graph
    return (ghash(nauty_graph), permutation[1:nv(cluster)])

end


function symmetric_pruning(cluster_prehash::AbstractNLCECluster)
    (
        hash(
            sort(([
                translational_form(
                                   cluster(underlying_lattice(cluster_prehash), Vector{Integer}(perm[vertices(cluster_prehash)])),
                ) for perm in permutations(underlying_lattice(cluster_prehash))
            ])),
        ),
        nothing,
    )
end
