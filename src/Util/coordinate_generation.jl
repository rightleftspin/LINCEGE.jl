"""
Generate every integer point in a cube of TWICE the given side length, centered around
the origin twice the side length is necessary to make sure a point exists at the origin
"""
function generate_cartesian_coordinates(dimension::Int, half_side_length::Int)
        # Forces the lattice to have a strict center point
        diameter = 2 * half_side_length + 1
        # Total number of coordinates for the entire lattice
        max_coords = diameter^dimension
        coords = repeat(0:(max_coords-1), 1, dimension)

        for dim = 0:(dimension-1)
                coords[:, dim+1] =
                        div.(coords[:, dim+1], diameter^dim) .% diameter .- half_side_length
        end

        coords
end

"""
Generates every point from the primitive lattice vectors, in a cube of
side length = 2 * maximum_order + 1 centered at the origin
"""
function generate_primitive_coordinates(primitive_vectors::AbstractMatrix{<:Real}, max_order::Int)

        unrotated_coords = generate_cartesian_coordinates(size(primitive_vectors, 1), max_order)
        # rotate and stretch standard cartesian cube coordinates into primitive lattice
        (unrotated_coords * stack(primitive_vectors), unrotated_coords)

end

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of real space coordinates
"""
function add_basis_coords(basis::AbstractMatrix{<:Real}, primitive_coordinates::AbstractMatrix{<:Real})
        n_basis = size(basis, 1)
        n_coords = size(primitive_coordinates, 1)
        dim = size(primitive_coordinates, 2)

        coords = Matrix{eltype(primitive_coordinates)}(undef, n_basis * n_coords, dim)

        for (i, shift) in enumerate(eachrow(basis))
                start_idx = (i - 1) * n_coords + 1
                end_idx = i * n_coords
                coords[start_idx:end_idx, :] .= primitive_coordinates .+ reshape(shift, 1, :)
        end

        return coords
end

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of sublattice coordinates
"""
function add_basis_sublattice(basis::AbstractMatrix{<:Real}, cartesian_coordinates::AbstractMatrix{<:Integer})
        n_basis = size(basis, 1)
        n_coords = size(cartesian_coordinates, 1)
        dim = size(cartesian_coordinates, 2)

        coords = Matrix{eltype(cartesian_coordinates)}(undef, n_basis * n_coords, dim + 1)

        for i in 1:n_basis
                start_idx = (i - 1) * n_coords + 1
                end_idx = i * n_coords
                coords[start_idx:end_idx, :] .= [cartesian_coordinates repeat([i], n_coords)]
        end

        return coords
end

"""
Generates the real space lattice from the expansion lattice, keeping track of which expansion lattice point 
is the origin of the given real space lattice point.
"""
function add_real_space_coords(
        tiling_units::AbstractVector{<:AbstractMatrix{<:Real}},
        expansion_coordinates::AbstractMatrix{<:Real},
        sublattice_coordinates::AbstractMatrix{<:Integer},
)
        n_coords = size(expansion_coordinates, 1)
        cs = []
        scs = []
        for (i, tiling_unit) in enumerate(tiling_units)
                start_idx = (i - 1) * n_coords + 1
                end_idx = i * n_coords
                n_basis_elems = size(tiling_unit, 1)
                push!(cs, add_basis_coords(tiling_unit, expansion_coordinates[start_idx:end_idx, :]))
                push!(scs, [add_basis_sublattice(tiling_unit, sublattice_coordinates[start_idx:end_idx, :]) repeat(1:n_coords, n_basis_elems)])
        end

        (hcat(cs...), hcat(scs...))
end

"Find the center of a set of coordinates"
function find_center(coordinates::AbstractMatrix{<:Real})
        sum(coordinates, dims=1) / size(coordinates, 1)
end

