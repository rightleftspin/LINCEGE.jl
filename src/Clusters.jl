
abstract type AbstractCluster end


"""
This is the cluster struct that powers the entire NLCE algorithm. It is designed to be very general so
it deals with many use cases.
"""
# Has underlying lattice, cluster expansion, vertex labeled, edge labeled
struct Cluster{U, C, V, E} <: AbstractCluster

    "Coordinates in the cluster, all index integers are in reference to this vector of coordinates"
    coordinates::AbstractVector{<:AbstractVector{<:Real}}
    "Start point of the cluster for growing subclusters, where elements of the vector point to specific coordinates in the vector of coordinates above"
    start::AbstractVector{<:Integer}
    "Bonds between sites where index in the outer vector is site 1 and index in the inner vector is site 2"
    adj_list::AbstractVector{<:AbstractVector{<:Integer}}
    "Bonds between sites where the first index is a choice of representation of the cluster, second index is site 1 and third index is site 2.
    The first two representation choices are reserved for isomorphic and translational hashing respectively.
    Since there are no self loops, the diagonal of the first adjacency matrix contains the labels of each site"
    adj_matrices::AbstractArray{<:Integer,3}

    "Underlying cluster for the given cluster, or nothing if it is the underlying lattice"
    underlying_cluster::Union{AbstractCluster, Nothing}
    "Vertices that are a part of the cluster in reference to an underlying lattice"
    vertices::Union{AbstractVector{<:Integer}, Nothing}

    # Below are the fields for the cluster expansion as opposed to the site expansion
    "Bundles of coordinates where each inner vector represents a single site in the super cluster in the case of a cluster expansion"
    coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing}
    "Bonds between sites in the super lattice. Here, the index integers are instead in reference to bundles in the coordinate bundles"
    super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing}

    function Cluster(
        coordinates::AbstractVector{<:AbstractVector{<:Real}},
        start::AbstractVector{<:Integer},
        adj_list::AbstractVector{<:AbstractVector{<:Integer}},
        adj_matrices::AbstractArray{<:Integer,3},
        underlying_cluster::Union{AbstractCluster, Nothing},
        vertices::Union{AbstractVector{<:Integer}, Nothing},
        coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
        super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
        has_underlying_cluster::Bool,
        cluster_expansion::Bool,
        vertex_labeled::Bool,
        edge_labeled::Bool,
    )
        if !has_underlying_cluster
            # Ensuring the cluster makes sense, but only at the
            # highest level so that you don't do a lot of redundant
            # checks
            n, _n = size(adj_matrices[1, :, :])
            m = length(adj_list)

            # Raising errors if it does not
            @assert n == _n "Adjacency Matrix needs to be square"
            @assert n == m "Adjacency Matrix and Adjacency list have different number of vertices"
        end

        return new{has_underlying_cluster, cluster_expansion, vertex_labeled, edge_labeled}(
            coordinates,
            start,
            adj_list,
            adj_matrices,
            underlying_cluster,
            vertices,
            coordinate_bundles,
            super_adj_list,
        )
    end
end

"""
Basic constructor for site expansion lattice with edge weights and vertex labels
"""
function Cluster(
    coordinates::AbstractVector{<:AbstractVector{<:Real}},
    start::AbstractVector{<:Integer},
    adj_matrices::AbstractArray{<:Integer,3},
    vertex_labeled::Bool,
    edge_labeled::Bool,
)

    Cluster(
        coordinates,
        start,
        adj_matrix_to_adj_list(adj_matrices[1, :, :]),
        adj_matrices,
        nothing,
        nothing,
        nothing,
        nothing,
        false,
        false,
        vertex_labeled,
        edge_labeled,
    )
end

# Constructor for cluster expansion lattice with edge weights and vertex labels

#TODO: Move the math heavy parts into util.jl
"""
Constructor for site expansion lattice, takes in a basis, primitive primitive vectors
and the maximum order and generates a cluster that represents the lattice.
"""
function Cluster(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer;
    basis_colors::AbstractVector{<:Integer} = repeat([1], length(basis)),
)

    coordinates, sublattice_coords, colors, start_points =
        generate_coordinates(basis, primitive_vectors, max_order, basis_colors)

    adj_matrices = zeros(Int, 2, length(coordinates), length(coordinates))
    directions::Vector{Vector{Real}} = []

    for (index_coord, coord) in enumerate(coordinates)
        adj_matrices[1, index_coord, index_coord] = colors[index_coord]
        for (index_distance, distance) in enumerate(neighborhood)
            equal_distance = n -> sqrt(sum((coord - n[2]) .^ 2)) ≈ distance
            # Find all neighbors equal to the current distance but after the last distance
            for (index_neighbor, neighbor) in
                filter(equal_distance, collect(enumerate(coordinates)))
                direction = neighbor - coord
                # Check the direction of the bond, ie, along which axis
                if findfirst(≈(direction), directions) != nothing
                    adj_matrices[1, index_coord, index_neighbor] = index_distance
                    adj_matrices[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                else
                    append!(directions, [direction])
                    adj_matrices[1, index_coord, index_neighbor] = index_distance
                    adj_matrices[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                end
            end
        end
    end

    Cluster(
        sublattice_coords,
        start_points,
        adj_matrices,
        (length(unique(basis_colors)) > 1),
        (length(neighborhood) > 1),
    )
end

"""
Takes the underlying cluster and returns a subcluster of it
"""
function Cluster(
    underlying_cluster::Cluster{<:Any, false, V, E},
    vertices::AbstractVector{<:Integer},
    ) where {V, E}
    Cluster(
        coordinates(underlying_cluster, vertices),
        Vector(1:length(vertices)),
        adj_matrix_to_adj_list(weighted_adjacency_matrix(underlying_cluster)[vertices, vertices]),
        adjacency_matrices(underlying_cluster)[:, vertices, vertices],
        underlying_cluster,
        vertices,
        nothing,
        nothing,
        true,
        false,
        V,
        E,
    )
end

begin # Standard Access functions
    nv(cluster::Cluster) = length(cluster.coordinates)
    coordinates(cluster::Cluster, vertices::Union{Integer, AbstractVector}) = cluster.coordinates[vertices]
    all_coordinates(cluster::Cluster) = cluster.coordinates
    start(cluster::Cluster) = cluster.start
    neighbors(cluster::Cluster{<:Any, false, <:Any, <:Any}, vertices::Union{Integer,AbstractVector}) = cluster.adj_list[vertices]
    adjacency_list(cluster::Cluster{<:Any, false, <:Any, <:Any}) = neighbors(cluster, 1:nv(cluster))
    edge_list(cluster::Cluster{<:Any, false, <:Any, <:Any}) = adj_matrix_to_edge_list(weighted_adjacency_matrix(cluster))
    adjacency_matrices(cluster::Cluster) = cluster.adj_matrices
    weighted_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[1, :, :] - diagm(labels(cluster))
    direction_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[2, :, :]
    label(cluster::Cluster, vertices::Union{Integer,AbstractVector}) = diag(weighted_adjacency_matrix(cluster))[vertices]
    labels(cluster::Cluster) = diag(adjacency_matrices(cluster)[1, :, :])
    underlying_cluster(cluster::Cluster{true, <:Any, <:Any, <:Any}) = cluster.underlying_cluster
    vertices(cluster::Cluster) = Vector(1:nv(cluster))
    underlying_vertices(cluster::Cluster{true, <:Any, <:Any, <:Any}) = cluster.vertices
    Base.show(io::IO, cluster::Cluster) = print(io, "Cluster with $(nv(cluster)) vertices\nCoordinates: $(all_coordinates(cluster))\nBonds: $(edge_list(cluster))")
end

begin # Hashing Functions
    """
    Takes a vertex colored cluster and finds the
    translationally invariant hash of it.

    Inputs:
          cluster: takes a cluster struct, uses the direction weights matrix in the cluster

    Output:
          Tuple of cluster hash and nothing
    """
    function translational_pruning(cluster::Cluster{true, false, true, <:Any})
        perm = sortperm(all_coordinates(cluster))
        form = vec(
            sum(
                weight -> 2^weight,
                direction_adjacency_matrix(cluster)[perm, :],
                dims = 2,
            )
        )
        (hash(form + (1 .// (labels(cluster)[perm] .+ 1))), nothing)
    end

    """
    Takes an uncolored cluster and finds the
    translationally invariant hash of it.

    Inputs:
          cluster: takes a cluster struct, uses the direction weights matrix in the cluster

    Output:
          Tuple of cluster hash and nothing
    """
    function translational_pruning(cluster::Cluster{true, false, false, <:Any})
        (hash(vec(
            sum(
                weight -> 2^weight,
                direction_adjacency_matrix(cluster)[sortperm(all_coordinates(cluster)), :],
                dims = 2,
            )
        )), nothing)
    end

    """
    Finds the canonical ordering of an edge-labeled, vertex colored cluster
    using Nauty, this function rearranges the cluster according to the permutation
    used by Nauty

    Inputs:
          cluster: takes in a cluster struct, uses the
                weighted_adjacency_matrix inside it

    Output:
          Tuple of cluster hash and permutation from nauty
    """
    function isomorphic_pruning(cluster::Cluster{true, false, <:Any, true})

        # This is to get edge-labels to work
        nauty_labels = vcat(labels(cluster),
                            zeros(Int64, fld(count(>(1), weighted_adjacency_matrix(cluster)), 2)))

        unweighted_adjacency_matrix = zeros(Int64, length(nauty_labels), length(nauty_labels))

        current_aux_vert = nv(cluster) + 1
        for j = 1:size(weighted_adjacency_matrix(cluster),2)
            for i = j:size(weighted_adjacency_matrix(cluster),1)

                if (weighted_adjacency_matrix(cluster)[i, j] == 1)
                    # Add an edge if there is already an edge
                    unweighted_adjacency_matrix[i, j] = 1
                    unweighted_adjacency_matrix[j, i] = 1
                elseif (weighted_adjacency_matrix(cluster)[i, j] !=0)
                    # Add an edge to the auxilary vertex here
                    unweighted_adjacency_matrix[i, current_aux_vert] = 1
                    unweighted_adjacency_matrix[j, current_aux_vert] = 1
                    unweighted_adjacency_matrix[current_aux_vert, i] = 1
                    unweighted_adjacency_matrix[current_aux_vert, j] = 1
                    # Color the aux vertex the same as the edge
                    nauty_labels[current_aux_vert] = weighted_adjacency_matrix(cluster)[i, j]
                    current_aux_vert += 1
                end
            end
        end

        nauty_graph =
            NautyGraph(unweighted_adjacency_matrix, nauty_labels)
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

    """
    Finds the canonical ordering of a vertex colored cluster
    using Nauty, this function rearranges the cluster according to the permutation
    used by Nauty

    Inputs:
          cluster: takes in a cluster struct, uses the
                weighted_adjacency_matrix inside it

    Output:
          Tuple of cluster hash and permutation from nauty
    """
    function isomorphic_pruning(cluster::Cluster{true, false, <:Any, false})

        nauty_graph =
            NautyGraph(weighted_adjacency_matrix(cluster), labels(cluster))
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

   # TODO: Need to fix the symmetric pruning function
   # function symmetric_pruning(cluster_prehash::AbstractNLCECluster)
   #     (
   #         hash(
   #             sort(([
   #                 translational_form(
   #                     cluster(underlying_lattice(cluster_prehash), Vector{Integer}(perm[vertices(cluster_prehash)])),
   #                 ) for perm in permutations(underlying_lattice(cluster_prehash))
   #                     ])),
   #         ),
   #         nothing,
   #     )
   # end

end

begin # Growing functions to get subclusters from a cluster
# TODO: Document all these functions, then optimize or parallelize them in some way

"""
This is step one of the NLCE pipeline. In this step, the algorithm takes in
a lattice and the order that clusters should be generated till.
Using this information, the algorithm recursively generates an array of
clusters that are all subclusters of specified order of the lattice.
"""
function grow(underlying_cluster::Cluster{<:Any, false, <:Any, <:Any}, max_order::Integer)
    out_array::Vector{AbstractCluster} = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in start(underlying_cluster)
        init_neighbors::Set{Int} = Set(
            collect(
                filter(neighbor -> !(neighbor in guarding_set), neighbors(underlying_cluster, vertex)),
            ),
        )
        vertices = [vertex]
        _grow_from_site(
            underlying_cluster,
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
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array. This adds all subclusters up
until the appropriate order, which is why it is only for the lattice
and not for individual clusters, which behave a bit differently

Inputs:
      lattice: cluster with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the order that
      the subclusters will be generated at

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the
      output is the out_array that gets added to.
"""
function _grow_from_site(
    lattice::Cluster{false, false, <:Any, <:Any},
    max_order::Integer,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::Set{V},
    guarding_set::Set{V},
    out_array::AbstractVector{<:AbstractCluster},
) where {V<:Integer}

    push!(out_array, Cluster(lattice, subcluster_vertices))

    if length(subcluster_vertices) == max_order
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(lattice, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            lattice,
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
        if (nv(lattice) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
Grows clusters from the given vertices of the underlying_cluster. Calls the code above, but
for each order until the max order in the cluster
"""
function grow(underlying_cluster::Cluster{true, <:Any, <:Any, <:Any})
    out_array::Vector{AbstractCluster} = Vector()

    for max_order = 1:(nv(underlying_cluster)-1)
        append!(out_array, grow(underlying_cluster, max_order))
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array.

Inputs:
      underlying_cluster: Graph with coordinates as vertex labels,
      vertex colors, and edge weights

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the
      output is the out_array that gets added to.
"""
function _grow_from_site(
    underlying_cluster::Cluster{true, false, <:Any, <:Any},
    max_order,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::AbstractSet{V},
    guarding_set::AbstractSet{V},
    out_array::AbstractVector{<:AbstractCluster},
) where {V<:Integer}

    if length(subcluster_vertices) == max_order
        push!(out_array, Cluster(underlying_cluster, subcluster_vertices))
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(underlying_cluster, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            underlying_cluster,
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
        if (nv(underlying_cluster) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

end

