# Need to do
#       grow:
#               - add type specifications for the inputs
#               - document the function itself
#               - convert input type for starting vertices to coordinates instead of Int64
#
#       grow_from_site:
#               - add type specifications for the inputs
#               - document the function itself
#               - optimize the function, it really needs it

"""
Main function in step one of the pipeline. Grows clusters from the given
vertices of the underlying_cluster.

Inputs: 
      underlying_cluster: MetaGraph with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the highest order that
      the subclusters will be generated till

      starting_vertices: array of vertex labels (coordinates)
      that are the coordinates that the clusters will be grown from

Output:
      Array of MetaGraph objects that are subclusters of the input
      graphs
"""
function grow(
        underlying_cluster, 
        max_order::Int64, 
        starting_vertices::Vector{Int64}
    )
    
    out_array = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in starting_vertices
        neighbors::Set{Int} = Set(collect(filter(neighbor -> !(neighbor in guarding_set), collect(outneighbors(underlying_cluster, vertex)))))
        vertices = [vertex]
        grow_from_site(
            underlying_cluster, 
            max_order, 
            vertices, 
            neighbors, 
            guarding_set,
            out_array
        )
        push!(guarding_set, vertex)
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array.

Inputs: 
      underlying_cluster: MetaGraph with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the highest order that
      the subclusters will be generated till

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the 
      output is the out_array that gets added to.
"""
function grow_from_site(
        underlying_cluster, 
        max_order::Int, 
        subcluster_vertices::Vector{Int}, 
        neighbors::Set{Int}, 
        guarding_set::Set{Int},
        out_array
    )
    
    if size(subcluster_vertices)[1] == max_order
        push!(out_array, induced_subgraph(underlying_cluster, subcluster_vertices)[1])
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(neighbors)
        neighbor = pop!(neighbors)
        append!(subcluster_vertices, neighbor) 

        new_neighbors = copy(neighbors)

        for vertex in outneighbors(underlying_cluster, neighbor)
            if (!(vertex in subcluster_vertices) & 
                !(vertex in new_guarding_set) & 
                !(vertex in new_neighbors))

                push!(new_neighbors, vertex)
            end
        end

        if grow_from_site(underlying_cluster, 
                          max_order, 
                          subcluster_vertices, 
                          new_neighbors, 
                          new_guarding_set,
                          out_array
                         )
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return(has_int_leaf)
        end
        push!(new_guarding_set, neighbor)
        if (nv(underlying_cluster) - length(new_guarding_set)[1]) < max_order
            return(has_int_leaf)
        end
    end
    return(has_int_leaf)
end
