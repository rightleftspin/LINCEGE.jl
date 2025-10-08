abstract type AbstractAdjacencyMatrix end

nv(g::AbstractAdjacencyMatrix) = size(g.adj_matrix, 1)
nw(g::AbstractAdjacencyMatrix) = length(g.weight_info)

labels(g::AbstractAdjacencyMatrix) = diag(g.adj_matrix)
adj_matrix(g::AbstractAdjacencyMatrix) = tril(g.adj_matrix, -1) + tril(g.adj_matrix, 1)

edge_list(g::AbstractAdjacencyMatrix) = adj_matrix_to_edge_list(adj_matrix(g))


Base.show(io::IO, g::AbstractAdjacencyMatrix) = print(
        io,
        "Graph with $(nv(g)) vertices, $(length(edge_list(g))) bonds, and $(nw(g)) unique weights",
)
