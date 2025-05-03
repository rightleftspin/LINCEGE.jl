"""
General utility functions for a variety of parts of the NLCE process, these
will generally be math heavy functions that are used often
"""

using LinearAlgebra

"""
Generates a lattice from the given basis and primitive lattice vectors, colors according to
the appropriate coloring
"""
function generate_coordinates(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    max_order::Integer,
    basis_colors::AbstractVector{<:Integer},
)
    primitive_lattice, cartesian_lattice = generate_primitive_lattice(primitive_vectors, max_order)
    coordinates = add_basis_coords(basis, primitive_lattice)
    sublattice_points = add_basis_sublattice(basis, cartesian_lattice)
    colors = repeat(basis_colors, size(primitive_lattice)[2])
    
    sorted_coord_perm = sortperm(coordinates)
    
    (coordinates[sorted_coord_perm],
     sublattice_points[sorted_coord_perm],
     colors[sorted_coord_perm],
     findfirst.(isapprox.(basis), (coordinates[sorted_coord_perm], ))
    )
end

"""
Generates a clustered lattice based on the given superlattice coordinates
"""
function generate_sub_coordinates(
    sup_coordinates::AbstractVector{<:AbstractVector{<:Real}},
    sup_sublattice_coordinates::AbstractVector{<:AbstractVector{<:Real}},
    sub_basis::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
    sub_basis_colors::AbstractVector{<:AbstractVector{<:Integer}},
)


    sub_coordinates::Vector{Vector{Real}} = []
    sub_colors::Vector{Int64} = []
    coordinate_bundles::Vector{Vector{Int64}} = []
    count = 1
    for (ind, coord) in enumerate(sup_sublattice_coordinates)
        append!(sub_coordinates, [sup_coordinates[ind] + elem for elem in sub_basis[coord[end]]])
        append!(sub_colors, sub_basis_colors[coord[end]])
        push!(coordinate_bundles, collect(count:(count + length(sub_basis[coord[end]]) - 1)))
        count += length(sub_basis[coord[end]])
    end

    (sub_coordinates,
     sub_colors,
     coordinate_bundles
     )
end

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of sublattice coordinates
"""
function add_basis_sublattice(basis, cartesian_lattice)
     collect(Iterators.flatten([[[site..., i] for i in 1:length(basis)] for site in eachcol(cartesian_lattice)]))
end

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of real space coordinates
"""
function add_basis_coords(basis, lattice)
    collect(Iterators.flatten([[site + elem for elem in basis] for site in eachcol(lattice)]))
end

"""
Generates every point from the primitive lattice vectors, in a cube of 
side length = 2 * maximum_order + 1 centered at the origin
"""
function generate_primitive_lattice(primitive_vectors, max_order)
   
    unrotated_coords = generate_cartesian_coordinates(length(primitive_vectors), max_order)
    # rotate and stretch standard cartesian cube coordinates into primitive lattice
    (stack(primitive_vectors) * unrotated_coords, unrotated_coords)

end

"""
Generate every integer point in a cube of TWICE the given side length, centered around 
the origin twice the side length is necessary to make sure a point exists at the origin
"""
function generate_cartesian_coordinates(dimension, half_side_length)
    # Forces the lattice to have a strict center point
    diameter = 2 * half_side_length + 1
    # Total number of coordinates for the entire lattice
    max_coords = diameter ^ dimension
    coords = repeat(transpose(0:max_coords - 1), dimension)
    
    for dim = 0:(dimension - 1)
        coords[dim + 1, :] = div.(coords[dim + 1, :], diameter ^ dim) .% diameter .- half_side_length
    end

    coords
end

#TODO: This could be much more efficient if I just constructed the adj list
#directly, do that eventually
"""
Creates just the adjacency list for a set of coordinates and neighborhoods
"""
function create_adj_list(coordinates, neighborhood)

    adj_matrix = zeros(Int, length(coordinates), length(coordinates))

    for (index_coord, coord) in enumerate(coordinates)
        for (index_distance, distance) in enumerate(neighborhood)
            equal_distance = n -> sqrt(sum((coord - n[2]) .^ 2)) ≈ distance
            # Find all neighbors equal to the current distance but after the last distance
            for (index_neighbor, neighbor) in
                filter(equal_distance, collect(enumerate(coordinates)))
                adj_matrix[index_coord, index_neighbor] = index_distance
            end
        end
    end

    adj_matrix_to_adj_list(adj_matrix)
end

"""
Create the adjacency matrices given a set of coordinates, colors and neighborhood.
This finds directions, and labels further nearest neighbor bonds in the
corresponding adjacency matrix.
"""
function create_adj_matrices(coordinates, colors, neighborhood)

    adj_matrices = zeros(Int, 2, length(coordinates), length(coordinates))
    directions::Vector{Vector{Real}} = []

    for (index_coord, coord) in enumerate(coordinates)
        adj_matrices[1, index_coord, index_coord] = colors[index_coord]
        for (index_distance, distance) in enumerate(neighborhood)
            equal_distance = n -> sqrt(sum((coord - n[2]) .^ 2)) ≈ distance
            # Find all neighbors equal to the current distance but after the last distance
            for (index_neighbor, neighbor) in
                filter(equal_distance, collect(enumerate(coordinates)))
                direction = neighbor - coord
                # Check the direction of the bond, ie, along which axis
                if findfirst(≈(direction), directions) != nothing
                    adj_matrices[1, index_coord, index_neighbor] = index_distance
                    adj_matrices[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                else
                    append!(directions, [direction])
                    adj_matrices[1, index_coord, index_neighbor] = index_distance
                    adj_matrices[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                end
            end
        end
    end

    adj_matrices
end

"""
Converts adjacency list for a graph to an adjacency matrix
"""
function adj_list_to_adj_matrix(
    adj_list::AbstractVector{<:AbstractVector{<:Integer}},
    values::AbstractVector{<:AbstractVector{<:Integer}},
)

    # Find the number of vertices in the graph
    number_vertices = length(adj_list)
    # Initialize the empty adjacency matrix
    adj_matrix::Matrix{Int64} = zeros(number_vertices, number_vertices)

    # Loop over all the vertices in the adjacency list
    for (vertex, neighbors) in enumerate(adj_list)
        # Loop over all their neighbors
        for (neighbor_index, neighbor) in enumerate(neighbors)
            # Add the corresponding value to the adjacency matrix
            adj_matrix[vertex, neighbor] = values[vertex][neighbor_index]
        end
    end

    adj_matrix
end

function adj_list_to_adj_matrix(adj_list::AbstractVector{<:AbstractVector{<:Integer}})

    # Wrapper function for an adjacency list without weights. This function
    # will instead add a weight of 1 for every edge
    values = deepcopy(adj_list)
    for i = 1:length(values)
        values[i] .= 1
    end

    adj_list_to_adj_matrix(adj_list, values)
end

"""
Converts adjacency matrix for a graph to an adjacency list
"""
function adj_matrix_to_adj_list(adj_matrix::AbstractMatrix{<:Integer})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]

    adj_list::Vector{Vector{Int64}} = []

    for i = 1:number_vertices
        neighbors = []
        for j = 1:number_vertices
            if i != j
                if adj_matrix[i, j] != 0
                    append!(neighbors, j)
                end
            end
        end
        push!(adj_list, neighbors)
    end

    adj_list

end

"""
Converts adjacency matrix for a graph to an edge list"""
function adj_matrix_to_edge_list(adj_matrix::AbstractMatrix{<:Integer})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]

    edge_list::Vector{Vector{Int64}} = []

    for i = 1:number_vertices
        for j = i:number_vertices
            if i != j
            if adj_matrix[i, j] != 0
                push!(edge_list, [i, j, adj_matrix[i, j]])
            end
            end
        end
    end

    edge_list

end

# TODO: This might need to be fixed? It seems right, I cant figure it out rn
function reindex_adj_list(underlying_adj_list::AbstractVector{<:AbstractVector{<:Integer}}, verts::AbstractVector{<:Integer})
    adj_matrix_to_adj_list(adj_list_to_adj_matrix(underlying_adj_list)[verts, verts])
end

# TODO: Fix this omg it is so slow
"""
Finds a permutation representation of a given group of transformations
on the given coordinates. 
    
"""
function find_permutations(coordinates, group, shifts)
    permutations::Vector{Vector{Union{Int, Nothing}}} = []

    for elem in group
        transformed_coords = [(elem * coord) for coord in coordinates]

#        perm_lengths = zeros(length(shifts))
#
#        for shift_i = 1:length(shifts)
#            perm_lengths[shift_i] = 
#
        max_perm = findfirst.(isapprox.([(coord + shifts[1]) for coord in transformed_coords], atol = 1e-5), (coordinates,))

        # This loop needs to be optimized, very parallelizeable
        for shift in shifts[2:end]
            perm = findfirst.(isapprox.([(coord + shift) for coord in transformed_coords], atol = 1e-5), (coordinates,))
            if length(unique(perm)) > length(unique(max_perm))
                max_perm = copy(perm)
            end
        end
        push!(
            permutations,
            max_perm
        )
    end

    println(length.(unique.(permutations)))
    
    permutations
end

"""
Resummation functions that take in thermodynamic properties
and return the resummed versions at various levels of correction.

NOTE: Please make sure to use at least two resummation techniques as
they will diverge at a further point than the bare NLCE summation.
This divergence between the resummation should be treated as the new
point to which the data is reasonable.
"""
