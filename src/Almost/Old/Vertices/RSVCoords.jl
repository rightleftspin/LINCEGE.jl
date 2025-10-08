
struct RSVCoords{C<:AbstractVector{<:AbstractVector{<:Real}},R<:AbstractVector{<:Int}}
        coordinates::C
        related_ev_coord::R
end

"Constructor for RSVCoords, generates all real space lattice coordinates from the expansion lattice coordinates"
function RSVCoords(tiling::Tiling, ev_coords::EVCoords)

        coordinates::Vector{Vector{Real}} = []
        related_ev_coord::AbstractVector{<:Int} = []
        for (ind, ev, sev) in zip(1:length(ev_coords), coordinates(ev_coords), sublattice_coordinates(ev_coords))
                append!(coordinates, shift_coords(tiling_units(tiling)[sev[end]], ev))
                append!(related_ev_coord, ind)
        end

        RSVCoords(
                coordinates,
                related_ev_coord
        )
end
