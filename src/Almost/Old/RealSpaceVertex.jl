struct RealSpaceVertex{C<:AbstractArray,R<:RealSpaceVertices,L<:Integer}
        coord::C
        label::L
        translational_label::L
end

function merge(v1::RealSpaceVertex, v2::RealSpaceVertex)

        @assert v1 == v2 "Cannot merge two RealSpaceVertex structs that do not refer to the same coordinate"
        @assert label(v1) == label(v2) "Labeling of the lattice violates translational symmetry, please choose a different labeling"

        RealSpaceVertex(
                coord(v1),
                label(v1),
                min(translational_label(v1), translational_label(v2)),
        )

end

function merge_coords(vs::AbstractSet{RealSpaceVertex}, v::RealSpaceVertex)
        # Should not encounter the exact same labeled coordinate twice,
        # truly bending the system with this definition of equality here
        if v in vs
                return push!(merge(pop!(vs, v), v))
        end

        push!(vs, v)
end

function reduce_coords(vs::AbstractVector{RealSpaceVertex})
        sort(collect(foldl(merge_coords, vs, init=Set{RealSpaceVertex}())))
end

label(v::RealSpaceVertex) = v.label
translation_label(v::RealSpaceVertex) = v.translational_label
