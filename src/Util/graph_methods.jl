function adj_mat_to_adj_list(adj_matrix::AbstractMatrix{<:Int})
        n_coords, _ = size(adj_matrix)
        adj_list = Vector{Vertices}(undef, n_coords)
        for i in 1:n_coords
                neighbors = Vertices()
                for j in 1:n_coords
                        if adj_matrix[i, j] > 0
                                neighbors = union(neighbors, Vertices(j))
                        end
                end
                adj_list[i] = neighbors
        end
        adj_list
end

function pairwise_distance_mat_to_adj_mat(pw_dist::AbstractMatrix{<:Real}, distances::AbstractVector{<:Real})
        n_coords, _ = size(pw_dist)
        index_matrix = Matrix{Int}(undef, size(pw_dist))
        dist_map = Dict(d => i for (i, d) in enumerate(distances))
        for i in 1:n_coords
                for j in 1:n_coords
                        index_matrix[i, j] = get(dist_map, pw_dist[i, j], 0)
                end
        end

        index_matrix
end
