struct EVCoords{C<:AbstractVector{<:AbstractVector{<:Real}},SC<:AbstractVector{<:AbstractVector{<:Int}}}
        coordinates::C
        sublattice_coordinates::SC
end

"Constructor for EVCoords, generates all expansion lattice coordinates from a given tiling"
function EVCoords(tiling::Tiling, max_order::Int)

        centers = find_centers(tiling)

        primitive_lattice, unrotated_primitive_lattice = generate_primitive_lattice(primitive_vectors(tiling), max_order)

        EVCoords(
                add_basis_coords(centers, primitive_lattice),
                add_basis_sublattice(centers, unrotated_primitive_lattice)
        )

end

Base.length(ev_coords::EVCoords) = length(ev_coords.coordinates)
Base.getindex(ev_coords::EVCoords, i::Int) = ev_coords.coordinates[i]
Base.getindex(ev_coords::EVCoords, i::Vector{Int}) = ev_coords.coordinates[i]
Base.getindex(ev_coords::EVCoords, i::ExpansionVertices) = ev_coords.coordinates[i]

sublattice_coordinate(ev_coords::EVCoords, i::Int) = ev_coords.sublattice_coordinates[i]
sublattice_coordinate(ev_coords::EVCoords, i::Vector{Int}) = ev_coords.sublattice_coordinates[i]
sublattice_coordinate(ev_coords::EVCoords, i::ExpansionVertices) = ev_coords.sublattice_coordinates[i]
