"""
This is the cluster struct that powers the entire NLCE algorithm. It is designed to be very general so
it deals with many use cases.
"""
#vertex labeled, edge labeled
struct Cluster{V,E}

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

        return new{vertex_labeled,edge_labeled}(adj_list, connections, adj_matrices)
    end
end

"""
Takes the underlying cluster and returns a subcluster of it.
"""
function Cluster(
    underlying_cluster::Cluster{V,E},
    vertices::AbstractVector{<:Integer},
) where {V,E}

    adj_list, connections, adj_matrices =
        reindex_to_subcluster(underlying_cluster, vertices)

    Cluster(adj_list, connections, adj_matrices, V, E)
end
