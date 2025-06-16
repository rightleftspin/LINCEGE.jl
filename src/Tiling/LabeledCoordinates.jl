struct RealSpaceVertex{C<:AbstractArray,L<:Integer,E<:ExpansionVertices}
        coord::C
        label::L
        translational_label::L
        expansion_vertices::E
end

function merge(l1::LabeledCoordinate, l2::LabeledCoordinate)

        # label for these two should be the same, if they are
        # related by translation, translational label just
        # devolves down to lowest possible label for all
        # coordinates
        LabeledCoordinate(
                l1.coord,
                l1.label,
                min(l1.translational_label, l2.translational_label),
                vcat(l1.expansion_vertex, l2.expansion_vertex)
        )

end

function merge_coords(ls::AbstractSet{LabeledCoordinate}, l::LabeledCoordinate)
        # Should not encounter the exact same labeled coordinate twice,
        # truly bending the system with this definition of equality here
        if l in ls
                return push!(merge(pop!(ls, l), l))
        end

        push!(ls, l)
end

function reduce_coords(ls::AbstractVector{LabeledCoordinate})
        sort(collect(foldl(merge_coords, ls, init=Set{LabeledCoordinate}())))
end


dimension(l::LabeledCoordinate) = length(l.coord)
coord(l::LabeledCoordinate) = l.coord
label(l::LabeledCoordinate) = l.label
translation_label(l::LabeledCoordinate) = l.translational_label
expansion_vertex(l::LabeledCoordinate) = l.expansion_vertex

Base.hash(l::LabeledCoordinate, h::UInt) = hash(l.coord, h)

Base.isequal(l1::LabeledCoordinate, l2::LabeledCoordinate) = (l1.coord == l2.coord)
