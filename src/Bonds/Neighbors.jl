struct Neighbors{A<:AbstractVector{<:Vertices}}
        adjacency_list::A
end

function Neighbors(tiling::Tiling, coordinates::Coordinates)
        adj_matrix = exp_neighbors(tiling, coordinates)

        Neighbors(
                adj_mat_to_adj_list(adj_matrix)
        )

end

neighbors(n::Neighbors, vs::Vertices) = n.adjacency_list[vs]
