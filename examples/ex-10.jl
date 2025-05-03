"""
Example showing a cluster expansion on a square lattice
"""
module ex10

using NLCE

#sup_basis = [[0, 0, 0]]
#sub_basis = [[[0, 0, 0], [0, 1 / 4, 1 / 4], [1 / 4, 0, 1 / 4], [1 / 4, 1 / 4, 0]]]
#
#sup_primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]
#sup_neighborhood = [sqrt(2) / 2]
#sub_neighborhood = [sqrt(2) / 4]

#sup_basis::Vector{Vector{Float64}} = [[0, 0]]
#sub_basis::Vector{Vector{Vector{Float64}}} = [[[0, 0], [1, 0], [0, 1], [1, 1]]]
#sup_primitive_vec::Vector{Vector{Float64}} = [[2, 0], [0, 2]]
#
#sup_neighborhood::Vector{Float64} = [2]
#sub_neighborhood::Vector{Float64} = [1]
#
#


max_order = 4
ssl_expansion_basis = [[0, 0], [2, 0]]
ssl_expansion_primitive_vectors = [[4, 0], [-2, 2]]
ssl_expansion_neighbors = [2]

ssl_struct_per_basis = [[[-0.5, 0], [0.5, 0]], [[0, 0.5], [0, -0.5]]]
ssl_colors = [[1, 2], [1, 2]]
ssl_neighbors = [1, sqrt(10) / 2]

println("start")

lattice = NLCE.Cluster(ssl_expansion_basis,
                       ssl_struct_per_basis,
                       ssl_expansion_primitive_vectors,
                       ssl_expansion_neighbors,
                       ssl_neighbors,
                       max_order)
println("lattice")

generated_clusters = NLCE.grow(lattice, max_order)
println("clusters")

iso_clusters = NLCE.prune(NLCE.isomorphic_pruning, Set(generated_clusters))
println("iso")

#for (hash, cluster) in iso_clusters
#    println(cluster[1])
#    println(cluster[2])
#end
#

propogated = NLCE.propogate(NLCE.isomorphic_pruning, iso_clusters)
println("propogated")

sums = NLCE.nlce_summation(propgated, max_order)
println(typeof(sums))


# TODO: Deal with single site multiplicity in here this could be put anywhere tbh


## Generating all the clusters using this information
#nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)
#
## Writing all the files to the corresponding folder, creating the folder
## if it does not exist
#filepath = "examples/outputs/ex-2/triangle_nn"
#mkpath(filepath)
#filename = filepath * "/triangle_nn"
#
## Write all the files in the "fortran" format
#write_to_file_fortran(nlce_clusters, filename, max_order)

end
