struct Lattice{EV<:EVCollection,RSV<:RSVCollection}

        expansion_vertices::EV
        real_space_vertices::RSV
end

function Lattice(tiling::Tiling, max_order::Int)

        ev_coords = EVCoords(tiling, max_order)

        rsv_coords = RSVCoords(tiling, ev_coords)

        real_space_vertices = RSVCollection(tiling, rsv_coords)

        expansion_vertices = EVCollection(tiling, ev_coords, real_space_vertices)

        Lattice(
                expansion_vertices,
                real_space_vertices
        )

end

function expansion_vertices(lattice::Lattice) end
