
function color_tiling_units(tiling_units::AbstractVector{<:AbstractMatrix{<:Real}}, primitive_vectors::AbstractMatrix{<:Real})

        pw_dir = pairwise_direction(vcat(tiling_units...))

        colors = find_colors(pw_dir, primitive_vectors)

        partition_colors(colors, tiling_units)
end

function pairwise_direction(coordinates::AbstractMatrix{<:Real})
        n_coords, dim = size(coordinates)
        pw_dir = Array{Real}(undef, n_coords, n_coords, dim)
        for i in 1:n_coords
                for j in 1:n_coords
                        pw_dir[i, j, :] = coordinates[j, :] .- coordinates[i, :]
                end
        end

        pw_dir
end

function find_colors(pw_dir::AbstractArray{<:Real,3}, primitive_vectors::AbstractMatrix{<:Real})
        n_coords, _, _ = size(pw_dir)
        colors = Vector(1:n_coords)
        prim_vector_math = transpose(primitive_vectors)

        for i in 1:n_coords
                for j in i:n_coords
                        direction = pw_dir[i, j, :]
                        int_coeffs = prim_vector_math \ direction
                        if all(isapprox.(int_coeffs, round.(int_coeffs)))
                                colors[i] = min(colors[i], colors[j])
                                colors[j] = min(colors[i], colors[j])
                        end

                end
        end

        colors
end

function partition_colors(colors::AbstractVector{<:Int}, tiling_units::AbstractVector{<:AbstractMatrix{<:Real}})
        partition_sizes = size.(tiling_units, 1)
        partitioned_colors = Vector{Vector{Int}}()
        start_idx = 1
        for partition_size in partition_sizes
                push!(partitioned_colors, colors[start_idx:partition_size+start_idx-1])
                start_idx += partition_size
        end

        partitioned_colors
end
