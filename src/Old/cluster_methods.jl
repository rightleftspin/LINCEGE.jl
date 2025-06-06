begin # Access functions, should take no computational effort
    nv(cluster::Cluster) = size(cluster.adj_matrices)[2]
    nsv(cluster::Cluster) = length(cluster.connections)

    vertex_labeled(cluster::Cluster{V,<:Any}) where {V} = V
    edge_labeled(cluster::Cluster{<:Any,E}) where {E} = E

    # Vertex Access Functions
    vertices(cluster::Cluster) = Vector(1:nv(cluster))
    super_vertices(cluster::Cluster) = Vector(1:nsv(cluster))
    vertex_label(cluster::Cluster, vertex::Integer) = adjacency_matrices[1, vertex, vertex]
    translation_label(cluster::Cluster, vertex::Integer) =
        adjacency_matrices[2, vertex, vertex]
    rev_connection(cluster::Cluster, vertex::Integer) =
        adjacency_matrices[3, vertex, vertex]

    all_vertex_labels(cluster::Cluster) = diag(adjacency_matrices(cluster)[1, :, :])
    all_translation_labels(cluster::Cluster) = diag(adjacency_matrices(cluster)[2, :, :])
    all_rev_connections(cluster::Cluster) = diag(adjacency_matrices(cluster)[3, :, :])

    # Edge Access Function
    adj_list(cluster::Cluster) = cluster.adj_list
    neighbors(cluster::Cluster, vertex::Integer) = adj_list(cluster)[vertex]
    bond_weight(cluster::Cluster, vertex1::Integer, vertex2::Integer) =
        weighted_adjacency_matrix(cluster)[vertex1, vertex2]
    bond_direction(cluster::Cluster, vertex1::Integer, vertex2::Integer) =
        direction_adjacency_matrix(cluster)[vertex1, vertex2]
    bond_sv(cluster::Cluster, vertex1::Integer, vertex2::Integer) =
        bond_sv_adjacency_matrix(cluster)[vertex1, vertex2]

    # Adjacency Matrix Access Functions
    adjacency_matrices(cluster::Cluster) = cluster.adj_matrices
    weighted_adjacency_matrix(cluster::Cluster) =
        adjacency_matrices(cluster)[1, :, :] - diagm(all_vertex_labels(cluster))
    direction_adjacency_matrix(cluster::Cluster) =
        adjacency_matrices(cluster)[2, :, :] - diagm(all_translation_labels(cluster))
    bond_sv_adjacency_matrix(cluster::Cluster) =
        adjacency_matrices(cluster)[3, :, :] - diagm(all_rev_connections(cluster))

    # Connection Access Functions
    connections(cluster::Cluster) = cluster.connections
    connection(cluster::Cluster, super_vertex::Integer) = connections(cluster)[super_vertex]

end

begin # Complex Access functions, might take some computational effort
    weighted_edge_list(cluster::Cluster) =
        adj_matrix_to_edge_list(weighted_adjacency_matrix(cluster))

    function reindex_to_subcluster(
        cluster::Cluster,
        super_vertices::AbstractVector{<:Integer},
    )
        # Sort them for translational hashing
        sorted_vertices = sort(unique(vcat(connections(cluster)[super_vertices]...)))

        new_adj_list = reindex_adj_list(adj_list(cluster), super_vertices)
        new_connections =
            reindex_connections(connections(cluster), super_vertices, sorted_vertices)
        new_adjacency_matrices = reindex_adjacency_matrices(
            adjacency_matrices(cluster),
            super_vertices,
            sorted_vertices,
        )

        new_adj_list, new_connections, new_adjacency_matrices
    end

    Base.show(io::IO, cluster::Cluster) = print(
        io,
        "Cluster with $(nv(cluster)) vertices and $(length(weighted_edge_list(cluster))) bonds. Super lattice contains $(nsv(cluster)) super vertices",
    )

    # Sets default hashing of a cluster to be the translationally invariant hash
    Base.hash(cluster::Cluster, h::UInt) = hash(translational_pruning(cluster), h)
    Base.isequal(cluster1::Cluster, cluster2::Cluster) =
        (translational_pruning(cluster1) == translational_pruning(cluster2))
end
