struct EVCollection{V<:EVCoords,A<:AbstractVector{<:ExpansionVertices},C<:AbstractVector{<:RealSpaceVertices}}
        ev_coords::V
        adjacency::A
        connections::C

end

function EVCollection(tiling::Tiling, ev_coords::EVCoords, real_space_vertex_collection::RSVCollection)

        EVCollection(
                ev_coords,
                adjacency,
                connections
        )
end

Base.getindex(evc::EVCollection, evs::ExpansionVertices) = evc.ev_coords[evs]
Base.sublattice_coordinates(evc::EVCollection, evs::ExpansionVertices) = sublattice_coordinate(evc.ev_coords, evs)
neighbors(evc::EVCollection, evs::ExpansionVertices) = union(evc.adjacency[evs][1], evc.adjacency[evs][2:end])
connections(evc::EVCollection, evs::ExpansionVertices) = union(evc.connections[evs][1], evc.connections[evs][2:end])
