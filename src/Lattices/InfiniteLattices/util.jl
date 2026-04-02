function generate_cartesian_coordinates(dimension::Int, half_side_length::Int)
        # Forces the lattice to have a strict center point
        r = -half_side_length:half_side_length
        hcat(vec(collect.(Iterators.product(fill(r, dimension)...)))...)
end

function generate_coordinates(max_order::Int, num_basis_elements::Int, dimension::Int)
        primitive_coordinates = generate_cartesian_coordinates(dimension, max_order)

        hcat([vcat(primitive_coordinates, repeat([i], size(primitive_coordinates, 2))') for i in 1:num_basis_elements]...)
end

function generate_coord_index(coordinates::Matrix{Int})
        coord_index = Dict{Vector{Int},Int}()
        for (i, col) in enumerate(eachcol(coordinates))
                coord_index[col] = i
        end
        coord_index
end

function generate_coord_index(coordinates::Matrix{Float64})
        coord_index = Dict{Vector{Float64},Int}()
        for (i, col) in enumerate(eachcol(coordinates))
                coord_index[round.(col, digits=6)] = i
        end
        coord_index
end

function generate_adj_matrix(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell)
        coord_index = generate_coord_index(coordinates)
        adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.bonds
                        if bond.site1 == col[end]
                                neighbor_coord = neighbor_site(bond, col)
                                ni = get(coord_index, neighbor_coord, nothing)
                                if ni !== nothing
                                        adj_matrix[index, ni] = bond.bond_type
                                        adj_matrix[ni, index] = bond.bond_type
                                end
                        end
                end
        end

        adj_matrix
end

function generate_adj_matrix(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell)
        coord_index = generate_coord_index(coordinates)
        adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.bonds
                        if bond.site1 == col[end-1:end]
                                neighbor_coord = neighbor_site(bond, col)
                                ni = get(coord_index, neighbor_coord, nothing)
                                if ni !== nothing
                                        adj_matrix[index, ni] = bond.bond_type
                                        adj_matrix[ni, index] = bond.bond_type
                                end
                        end
                end
        end

        adj_matrix
end

function generate_adj_matrix_weak(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell, unique_inds::Vector{Int})
        real_coords = shift_unit_cell(unit_cell, coordinates)[:, unique_inds]
        coord_index = generate_coord_index(real_coords)
        adj_matrix = zeros(Int, length(unique_inds), length(unique_inds))
        for coord in eachcol(coordinates)
                for bond in unit_cell.bonds
                        if bond.site1 == coord[end-1:end]
                                trans_ind = get(coord_index, round.(shift_unit_cell(unit_cell, coord), digits=6), nothing)
                                neighbor_coord = neighbor_site(bond, coord)
                                ni = get(coord_index, round.(shift_unit_cell(unit_cell, neighbor_coord), digits=6), nothing)
                                if ni !== nothing
                                        adj_matrix[trans_ind, ni] = bond.bond_type
                                        adj_matrix[ni, trans_ind] = bond.bond_type
                                end
                        end
                end
        end

        adj_matrix
end

function generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell)
        coord_index = generate_coord_index(coordinates)
        neighbor_list = fill(LatticeVertices{Int}(), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.bonds
                        if bond.site1 == col[end]
                                neighbor_index = get(coord_index, neighbor_site(bond, col), nothing)
                                if neighbor_index !== nothing
                                        neighbor_list[index] = union(neighbor_list[index], LatticeVertices(neighbor_index))
                                        neighbor_list[neighbor_index] = union(neighbor_list[neighbor_index], LatticeVertices(index))
                                end
                        end
                end
        end

        neighbor_list
end

function generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell)
        coord_index = generate_coord_index(coordinates)
        neighbor_list = fill(ExpansionVertices{Int}(), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.expansion_bonds
                        if bond.site1 == col[end]
                                neighbor_index = get(coord_index, neighbor_site(bond, col), nothing)
                                if neighbor_index !== nothing
                                        neighbor_list[index] = union(neighbor_list[index], ExpansionVertices(neighbor_index))
                                        neighbor_list[neighbor_index] = union(neighbor_list[neighbor_index], ExpansionVertices(index))
                                end
                        end
                end
        end

        neighbor_list
end

find_centers(coordinates::AbstractMatrix{Int}) = findall(col -> all(==(0), col[1:end-1]), eachcol(coordinates))

function generate_strong_connections(expansion_coordinates::AbstractMatrix{Int}, lattice_coordinates::AbstractMatrix{Int})

        connections_vec = fill(LatticeVertices{Int}(), size(expansion_coordinates, 2))
        for (i, coord) in enumerate(eachcol(expansion_coordinates))
                connection = findall(==(coord[1:end-1]), eachcol(lattice_coordinates[1:end-2, :]))
                connections_vec[i] = LatticeVertices(collect(connection))
        end
        connections_vec
end

# Source - https://stackoverflow.com/a/50900113
# Posted by Bogumił Kamiński
# Retrieved 2026-04-02, License - CC BY-SA 4.0
# Modified by adding a round function to deal
# with floating point numbers
function uniqueidx(x::AbstractArray{T}) where {T}
        uniqueset = Set{T}()
        ex = eachindex(x)
        idxs = Vector{eltype(ex)}()
        for i in ex
                xi = round.(x[i], digits=6)
                if !(xi in uniqueset)
                        push!(idxs, i)
                        push!(uniqueset, xi)
                end
        end
        idxs
end

function generate_weak_connections(expansion_coordinates::AbstractMatrix{Int}, lattice_coords::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell)
        connections_vec = fill(LatticeVertices{Int}(), size(expansion_coordinates, 2))

        real_coords = shift_unit_cell(unit_cell, lattice_coords)
        unique_inds = uniqueidx(Vector{Vector{Float64}}(eachcol(real_coords)))

        real_coords = real_coords[:, unique_inds]
        reduced_lattice_coordinates = lattice_coords[:, unique_inds]
        reverse_connections_vec = fill(ExpansionVertices{Int}(), size(reduced_lattice_coordinates, 2))

        for (i, coord) in enumerate(eachcol(expansion_coordinates))
                connection = Int[]
                for lcoord in eachcol(lattice_coordinates(unit_cell, coord))
                        con = findfirst(isapprox(lcoord), eachcol(real_coords))
                        if !isnothing(con)
                                push!(connection, con)
                        end
                end
                connections_vec[i] = LatticeVertices(connection)
                for lv in connection
                        reverse_connections_vec[lv] = union(reverse_connections_vec[lv], ExpansionVertices(i))
                end
        end

        masking_matrix = zeros(Int, size(reduced_lattice_coordinates, 2), size(reduced_lattice_coordinates, 2))
        for (i, evs1) in enumerate(reverse_connections_vec)
                for (j, evs2) in enumerate(reverse_connections_vec)
                        shared_evs = intersect(evs1, evs2)
                        if length(shared_evs) == 1
                                masking_matrix[i, j] = collect(shared_evs)[1]
                        end
                end
        end

        return connections_vec, reverse_connections_vec, masking_matrix, unique_inds
end




