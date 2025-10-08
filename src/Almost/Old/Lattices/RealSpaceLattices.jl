struct RealSpaceLattice{
        C<:AbstractVector{RealSpaceCoordinates},
}
        coordinates::C
end

function RealSpaceLattice(nlce_tiling::Tiling, max_order::Int)

end

Base.getindex(lattice::RealSpaceLattice, ind::Int) = lattice.coordinates[ind]

pairwise_distance(lattice::RealSpaceLattice) = pairwise(euclidean, vcat(coords(lattice)), dims=1)
coords(lattice::RealSpaceLattice) = coord.(lattice.coordinates)
labels(lattice::RealSpaceLattice) = label.(lattice.coordinates)
translational_labels(lattice::RealSpaceLattice) = translational_label.(lattice.coordinates)

exp_vertices(lattice::RealSpaceLattice, ind::Int) = exp_vertices(lattice[ind])

function matching_exp_vertex(lattice::RealSpaceLattice, i::Int, j::Int)
        # Sometimes, there are dangling bonds in the lattice, however, these are usually
        # far out enough to where you will never expand into then, so they can safely
        # be set to 0 and not be dealt with
        v = intersect(exp_vertices(lattice, i), exp_vertices(lattice, j))
        if isempty(v)
                return zero(eltype(v))
        end

        pop!(v)
end
