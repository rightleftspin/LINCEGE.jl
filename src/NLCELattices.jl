# Need to do
#       NLCELattice:
#                       - Raise errors if adj_list and adj_matrix don't make sense
#                       for their weights
#                       - Raise errors if adj_list_weights contains weights that are zeros
#                       - Raise errors if vertex_labels contain labels that are zero

"""
This is the lattice struct that powers the entire NLCE algorithm. This struct is an immutable 
"""

abstract type AbstractNLCELattice end

struct NLCELattice{D, W, L} <: AbstractNLCELattice
    number_vertices::Int
    vertex_labels::AbstractVector{<:Integer}, 
    adj_list::AbstractVector{<:AbstractVector{<:Integer}},
    adj_matrix::AbstractMatrix{<:Integer}, 
    adj_list_weights::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}}, 
    adj_matrix_weights::AbstractVector{<:AbstractMatrix{<:Integer}}

    # Basic Constructor
    function NLCELattice{D, W, L}(
                                  vertex_labels::AbstractVector{<:Integer}, 
                                  adj_list::AbstractVector{<:AbstractVector{<:Integer}},
                                  adj_matrix::AbstractMatrix{<:Integer}, 
                                  adj_list_weights::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}}, 
                                  adj_matrix_weights::AbstractVector{<:AbstractMatrix{<:Integer}},
                                  directed::Bool,
                                  edge_weighted::Bool,
                                  vertex_labeled::Bool
                                 )
        # Ensuring the lattice makes sense
        n, _n = size(adj_matrix)
        m = length(adj_matrix)
        lab = length(vertex_labels)

        # Raising errors if it does not
        @assert n == _n "Adjacency Matrix needs to be square"
        @assert n == m "Adjacency Matrix and Adjacency list have different number of vertices"
        @assert lab == m "Wrong number of vertex labels"

        return new{directed, edge_weighted, vertex_labeled}(n, vertex_labels, adj_list, adj_matrix, adj_list_weights, adj_matrix_weights)
    end

end

    # Adjacency List Constructor
    function NLCELattice(
                            vertex_labels::AbstractVector{<:Integer}, 
                            adj_list::AbstractVector{<:AbstractVector{<:Integer}},
                            adj_list_weights::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}}, 
                            directed::Bool,
                            edge_weighted::Bool,
                            vertex_labeled::Bool
                        )
        # Generate adj_matrix
        adj_matrix = adj_list_to_adj_matrix(adj_list)
        # Generate all weights for adj_matrix
        adj_matrix_weights = [adj_list_to_adj_matrix(adj_list, weights) for weights in adj_list_weights]

        return NLCELattice(vertex_labels, adj_list, adj_matrix, adj_list_weights, adj_matrix_weights, directed, edge_weighted, vertex_labeled)
    end

    # Adjacency Matrix Constructor
    function NLCELattice(
                            vertex_labels::AbstractVector{<:Integer}, 
                            adj_matrix::AbstractMatrix{<:Integer}, 
                            adj_matrix_weights::AbstractVector{<:AbstractMatrix{<:Integer}},
                            directed::Bool,
                            edge_weighted::Bool,
                            vertex_labeled::Bool
                        )
        # Generate adj_list
        adj_list = adj_matrix_to_adj_list(adj_matrix)
        # Generate all weights for adj_matrix
        adj_list_weights = [adj_matrix_to_adj_list(adj_matrix, weights) for weights in adj_matrix_weights]

        return NLCELattice(vertex_labels, adj_list, adj_matrix, adj_list_weights, adj_matrix_weights, directed, edge_weighted, vertex_labeled)
    end

begin #Required functions for grow part of the pipeline
    nv(lattice::NLCELattice) = lattice.number_vertices
    neighbors(lattice::NLCELattice, vertex::Int64) = @view lattice.adj_list[vertex]
    
    function cluster(lattice::NLCELattice, vertices::AbstractVector{<:Integer})


    end

end


