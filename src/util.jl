# Need to do
#       adj_list_to_adj_matrix:
#                       - optimize

"""
General utility functions for a variety of parts of the NLCE process, these
will generally be math heavy functions that are used often
"""

"""
Converts adjacency list for a graph to an adjacency matrix
"""
function adj_list_to_adj_matrix(adj_list::AbstractVector{<:AbstractVector{<:Integer}}, values::AbstractVector{<:AbstractVector{<:Integer}})

    # Find the number of vertices in the graph
    number_vertices = length(adj_list)
    # Initialize the empty adjacency matrix
    adj_matrix::Matrix{Int64} = zeros(number_vertices, number_vertices)

    # Loop over all the vertices in the adjacency list
    for (vertex, neighbors) in enumerate(adj_list)
        # Loop over all their neighbors
        for (neighbor_index, neighbor) in enumerate(neighbors)
            # Add the corresponding value to the adjacency matrix
            adj_matrix[vertex, neighbor] = values[vertex][neighbor_index]
        end
    end

    adj_matrix
end

function adj_list_to_adj_matrix(adj_list::AbstractVector{<:AbstractVector{<:Integer}})

    # Wrapper function for an adjacency list without weights. This function
    # will instead add a weight of 1 for every edge
    values = deepcopy(adj_list)
    for i in 1:length(values)
        values[i] .= 1
    end

    adj_list_to_adj_matrix(adj_list, values)
end

"""
Converts adjacency matrix for a graph to an adjacency list
"""
function adj_matrix_to_adj_list(adj_matrix::AbstractMatrix{<:Integer})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]
    
    adj_list::Vector{Vector{Int64}} = []
    
    for i in 1:number_vertices
        neighbors = []
        for j in 1:number_vertices
            if adj_matrix[i, j] != 0
                append!(neighbors, j)
            end
        end
        push!(adj_list, neighbors)
    end

    adj_list

end

function adj_matrix_to_adj_list(adj_matrix::AbstractMatrix{<:Integer}, adj_matrix_weights::AbstractMatrix{<:Integer})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]
    
    adj_list::Vector{Vector{Int64}} = []
    
    for i in 1:number_vertices
        neighbors = []
        for j in 1:number_vertices
            if adj_matrix[i, j] != 0
                append!(neighbors, adj_matrix_weights[i, j])
            end
        end
        push!(adj_list, neighbors)
    end

    adj_list

end
