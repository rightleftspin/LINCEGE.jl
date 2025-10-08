struct Coordinates{C<:AbstractMatrix{<:Real},SC<:AbstractMatrix{<:Int}}
        coordinates::C
        sublattice_coordinates::SC
end

"Constructor for Expansion Vertex Coordinates, generates all expansion lattice coordinates from a given tiling"
function Coordinates(tiling::Tiling, max_order::Int)
        centers = find_centers(tiling)

        primitive_lattice, unrotated_primitive_lattice = generate_primitive_lattice(translation_vectors(tiling), max_order)

        Coordinates(
                add_basis_coords(centers, primitive_lattice),
                add_basis_sublattice(centers, unrotated_primitive_lattice)
        )

end

"Constructor for Real Space Coordinates, generates all real space lattice coordinates from the expansion lattice coordinates"
function Coordinates(tiling::Tiling, ev_coords::Coordinates)

        coordinates, sublattice_coordinates = add_real_space_coords(tiling_units(tiling), all_coordinates(ev_coords), all_sublattice_coordinates(ev_coords))

        Coordinates(
                coordinates,
                sublattice_coordinates
        )
end

Base.length(coordinates::Coordinates) = size(coordinates.coordinates, 1)

coordinates(coordinates::Coordinates, v::Vertices) = coordinates.coordinates[v, :]
all_coordinates(coordinates::Coordinates) = coordinates.coordinates
sublattice_coordinates(coordinates::Coordinates, v::Vertices) = coordinates.sublattice_coordinates[v, :]
all_sublattice_coordinates(coordinates::Coordinates) = coordinates.sublattice_coordinates
pairwise_distance(coordinates::Coordinates) = pairwise(euclidean, all_coordinates(coordinates), dims=1)

