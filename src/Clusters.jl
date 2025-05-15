using NautyGraphs

"""
This is the cluster struct that powers the entire NLCE algorithm. It is designed to be very general so
it deals with many use cases.
"""
#vertex labeled, edge labeled
struct Cluster{V, E}

    "Bonds between sites in the super lattice. Here, the index integers are instead in reference to bundles in the coordinate bundles"
    adj_list::AbstractVector{<:AbstractVector{<:Integer}}

    "Connection between the sites in the adjacency list and the sites in the adjacency matrix, a 1 to 1 correspondence in the site expansion
    but becomes more complex in a cluster expansion"
    connections::AbstractVector{<:AbstractVector{<:Integer}}

    "Bonds between sites where the first index is a choice of representation of the cluster, second index is site 1 and third index is site 2.
    The first two representation choices are reserved for isomorphic and translational hashing respectively. The third rep. choice is reserved
    for connecting bonds to sites in the super lattice.
    Since there are no self loops, the diagonal of the first adjacency matrix contains the labels of each site, the diagonals of the
    second adjacency matrix contains the labels for each site that is used for translational invariance"
    adj_matrices::AbstractArray{<:Integer,3}

    function Cluster(
        adj_list::AbstractVector{<:AbstractVector{<:Integer}},
        connections::AbstractVector{<:AbstractVector{<:Integer}},
        adj_matrices::AbstractArray{<:Integer,3},
        vertex_labeled::Bool,
        edge_labeled::Bool,
    )

        # TODO: Put checks here to make sure the cluster makes sense

        return new{vertex_labeled, edge_labeled}(
            adj_list,
            connections,
            adj_matrices,
        )
    end
end

"""
Takes the underlying cluster and returns a subcluster of it.
"""
function Cluster(
    underlying_cluster::Cluster{V, E},
    vertices::AbstractVector{<:Integer},
    ) where {V, E}

    adj_list, connections, adj_matrices = reindex_to_subcluster(underlying_cluster, vertices)

    Cluster(
        adj_list,
        connections,
        adj_matrices,
        V,
        E,
    )
end

begin # Access functions, should take no computational effort
    nv(cluster::Cluster) = size(cluster.adj_matrices)[2]
    nsv(cluster::Cluster) = length(cluster.connections)

    vertex_labeled(cluster::Cluster{V, <:Any}) where {V} = V
    edge_labeled(cluster::Cluster{<:Any, E}) where {E} = E

    # Vertex Access Functions
    vertices(cluster::Cluster) = Vector(1:nv(cluster))
    super_vertices(cluster::Cluster) = Vector(1:nsv(cluster))
    vertex_label(cluster::Cluster, vertex::Integer) = adjacency_matrices[1, vertex, vertex]
    translation_label(cluster::Cluster, vertex::Integer) = adjacency_matrices[2, vertex, vertex]
    rev_connection(cluster::Cluster, vertex::Integer) = adjacency_matrices[3, vertex, vertex]

    all_vertex_labels(cluster::Cluster) = diag(adjacency_matrices(cluster)[1, :, :])
    all_translation_labels(cluster::Cluster) = diag(adjacency_matrices(cluster)[2, :, :])
    all_rev_connections(cluster::Cluster) = diag(adjacency_matrices(cluster)[3, :, :])

    # Edge Access Function
    adj_list(cluster::Cluster) = cluster.adj_list
    neighbors(cluster::Cluster, vertex::Integer) = adj_list(cluster)[vertex]
    bond_weight(cluster::Cluster, vertex1::Integer, vertex2::Integer) = weighted_adjacency_matrix(cluster)[vertex1, vertex2]
    bond_direction(cluster::Cluster, vertex1::Integer, vertex2::Integer) = direction_adjacency_matrix(cluster)[vertex1, vertex2]
    bond_sv(cluster::Cluster, vertex1::Integer, vertex2::Integer) = bond_sv_adjacency_matrix(cluster)[vertex1, vertex2]

    # Adjacency Matrix Access Functions
    adjacency_matrices(cluster::Cluster) = cluster.adj_matrices
    weighted_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[1, :, :] - diagm(all_vertex_labels(cluster))
    direction_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[2, :, :] - diagm(all_translation_labels(cluster))
    bond_sv_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[3, :, :] - diagm(all_rev_connections(cluster))

    # Connection Access Functions
    connections(cluster::Cluster) = cluster.connections
    connection(cluster::Cluster, super_vertex::Integer) = connections(cluster)[super_vertex]

end

begin # Complex Access functions, might take some computational effort
    weighted_edge_list(cluster::Cluster) = adj_matrix_to_edge_list(weighted_adjacency_matrix(cluster))

    function reindex_to_subcluster(cluster::Cluster, super_vertices::AbstractVector{<:Integer})
        # Sort them for translational hashing
        sorted_vertices = sort(unique(vcat(connections(cluster)[super_vertices]...)))

        new_adj_list = reindex_adj_list(adj_list(cluster), super_vertices)
        new_connections = reindex_connections(connections(cluster), super_vertices, sorted_vertices)
        new_adjacency_matrices = reindex_adjacency_matrices(adjacency_matrices(cluster), super_vertices, sorted_vertices)

        new_adj_list, new_connections, new_adjacency_matrices
    end

    Base.show(io::IO, cluster::Cluster) = print(io, "Cluster with $(nv(cluster)) vertices and $(length(weighted_edge_list(cluster))) bonds. Super lattice contains $(nsv(cluster)) super vertices")

    # Sets default hashing of a cluster to be the translationally invariant hash
    Base.hash(cluster::Cluster, h::UInt) = hash(translational_pruning(cluster), h)
    Base.isequal(cluster1::Cluster, cluster2::Cluster) = (translational_pruning(cluster1) == translational_pruning(cluster2))
end

begin # Hashing Functions

    """
    Takes a vertex labeled cluster and finds the translationally invariant hash of it.

    Please note that this function requires all the vertices to be sorted in dimensional order,
    ie. x, y, z, ...
    """
    # TODO: Add support for tiling here
    function translational_pruning(cluster::Cluster{true, <:Any})
        form = vec(
            sum(
                weight -> 2^weight,
                direction_adjacency_matrix(cluster),
                dims = 2,
            )
        )
        # TODO: Document this deviation from MVI
        (hash(form + (1 .// (all_vertex_labels(cluster) .+ 1))), nothing)
    end

    """
    Takes a vertex unlabeled cluster and finds the translationally invariant hash of it.

    Please note that this function requires all the vertices to be sorted in dimensional order,
    ie. x, y, z, ...
    """
    function translational_pruning(cluster::Cluster{false, <:Any})
        (hash(vec(
            sum(
                weight -> 2^weight,
                direction_adjacency_matrix(cluster),
                dims = 2,
            )
        ) + (1 .// (all_translation_labels(cluster) .+ 1))), nothing)
    end

    """
    Finds the canonical ordering of an edge labeled cluster using Nauty, this function
    returns the permutation to rearrange the cluster used by Nauty
    """
    function isomorphic_pruning(cluster::Cluster{<:Any, true})

        weighted_adj_mat = weighted_adjacency_matrix(cluster)
        # This is to get edge-labels to work
        nauty_labels = vcat(all_vertex_labels(cluster),
                            zeros(Int64, fld(count(>(1), weighted_adj_mat), 2)))

        unweighted_adjacency_matrix = zeros(Int64, length(nauty_labels), length(nauty_labels))

        current_aux_vert = nv(cluster) + 1
        for j = 1:size(weighted_adj_mat, 2)
            for i = j:size(weighted_adj_mat, 1)

                if (weighted_adj_mat[i, j] == 1)
                    # Add an edge if there is already an edge
                    unweighted_adjacency_matrix[i, j] = 1
                    unweighted_adjacency_matrix[j, i] = 1
                elseif (weighted_adj_mat[i, j] !=0)
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

        nauty_graph =
            NautyGraph(unweighted_adjacency_matrix, nauty_labels)
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
    function isomorphic_pruning(cluster::Cluster{<:Any, false})

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

begin # Growing functions to get subclusters from a cluster

"""
The algorithm takes in a cluster and the order that subclusters should be generated till.
Using this information, the algorithm recursively generates an array of vertices that
are all subclusters of specified order or lower of the cluster.
"""
function grow_lower(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

    out_array::Vector{Vector{Int}} = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in start
        init_neighbors::Set{Int} = Set(
            collect(
                filter(neighbor -> !(neighbor in guarding_set), neighbors(cluster, vertex)),
            ),
        )
        vertices = [vertex]
        _grow_from_site_lower(
            cluster,
            max_order,
            vertices,
            init_neighbors,
            guarding_set,
            out_array,
        )
        push!(guarding_set, vertex)
    end

    out_array
end


"""
This is step one of the NLCE pipeline. The algorithm takes in
a lattice and the order that clusters should be generated till.
Using this information, the algorithm recursively generates an array of
clusters that are all subclusters of specified order or lower of the lattice.
"""
function grow_exact(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

    out_array::Vector{Vector{Int}} = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in start
        init_neighbors::Set{Int} = Set(
            collect(
                filter(neighbor -> !(neighbor in guarding_set), neighbors(cluster, vertex)),
            ),
        )
        vertices = [vertex]
        _grow_from_site_exact(
            cluster,
            max_order,
            vertices,
            init_neighbors,
            guarding_set,
            out_array,
        )
        push!(guarding_set, vertex)
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order and outputs them into
the out_array. This adds all subclusters at or below the given max_order.
"""
function _grow_from_site_lower(
    cluster::Cluster,
    max_order::Integer,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::Set{V},
    guarding_set::Set{V},
    out_array::AbstractVector{<:AbstractVector{V}},
) where {V<:Integer}

    push!(out_array, deepcopy(subcluster_vertices))

    if length(subcluster_vertices) == max_order
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(cluster, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site_lower(
            cluster,
            max_order,
            subcluster_vertices,
            new_neighbors,
            new_guarding_set,
            out_array,
        )
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return (has_int_leaf)
        end
        push!(new_guarding_set, neighbor)

        if (nsv(cluster) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
Grows the subclusters from a specific site, up till specific order and outputs them into
the out_array. This adds all subclusters at the given max_order.
"""
function _grow_from_site_exact(
    cluster::Cluster,
    max_order::Integer,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::Set{V},
    guarding_set::Set{V},
    out_array::AbstractVector{<:AbstractVector{V}},
) where {V<:Integer}


    if length(subcluster_vertices) == max_order
        push!(out_array, deepcopy(subcluster_vertices))
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(cluster, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site_exact(
            cluster,
            max_order,
            subcluster_vertices,
            new_neighbors,
            new_guarding_set,
            out_array,
        )
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return (has_int_leaf)
        end
        push!(new_guarding_set, neighbor)

        if (nsv(cluster) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

end

