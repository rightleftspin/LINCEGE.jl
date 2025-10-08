# Neighbor functions take in Coordinates and return weighted adjacency matrices
# Label functions take in Coordinates and return vector of labels
struct Tiling{
        B<:AbstractVector{<:AbstractMatrix{<:Real}},
        T<:AbstractMatrix{<:Real},
        F,
}
        tiling_units::B
        translation_vectors::T
        neighbor_fn::F
        exp_neighbor_fn::F
        label_fn::F

        function Tiling(tiling_units::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}}, translation_vectors::AbstractVector{<:AbstractVector{<:Real}}, neighbor_fn, exp_neighbor_fn, label_fn)
                Tiling(
                        [collect(transpose(hcat(i...))) for i in tiling_units],
                        collect(transpose(hcat(translation_vectors...))),
                        neighbor_fn,
                        exp_neighbor_fn,
                        label_fn,
                )

        end
end

tiling_units(tiling::Tiling) = tiling.tiling_units
translation_vectors(tiling::Tiling) = tiling.translation_vectors
neighbors(tiling::Tiling, coordinates::Coordinates) = tiling.neighbor_fn(coordinates)
exp_neighbors(tiling::Tiling, coordinates::Coordinates) = tiling.exp_neighbor_fn(coordinates)
labels(tiling::Tiling, coordinates::Coordinates) = tiling.label_fn(coordinates)

"Find the center of each tiling unit in the tiling"
function find_centers(tiling::Tiling)
        find_center.(tiling_units(tiling))
end

# TODO: Create constructors for tiling in various ways. Base tiling is complete already
