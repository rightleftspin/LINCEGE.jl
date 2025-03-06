module ex8

using NLCE

basis = [[0, 0], [1, 0], [0, 1], [1, 1]]

colors = [1, 2, 2, 1]

primitive_vec = [[2, 0], [0, 2]]

neighborhood = [1]

max_order = 3

nlce_clusters = NLCE.site_color_NLCE(basis, colors, primitive_vec, neighborhood, max_order)

filepath = "examples/outputs/ex-8/ssl_dimer"
mkpath(filepath)
filename = filepath * "/ssl_dimer_$(max_order)"

NLCE.write_to_file_colors(nlce_clusters, filename)


# Rewriting clusters to be useful
for (cluster, multiplicities) in nlce_clusters
    intra_dimer = []
    for site in 0:(NLCE.nv(cluster) - 1)
        push!(intra_dimer, (2 * site, (2 * site + 1)))
    end

    inter_dimer = []
    #
    for bond in NLCE.edge_list(cluster)
        sorted_bond = sort(bond[[1, 2]]) .- 1
        site_1 = NLCE.all_coordinates(cluster)[sorted_bond[1] + 1]
        site_2 = NLCE.all_coordinates(cluster)[sorted_bond[2] + 1]

        sorted_sites = sort([site_1[[1, 2]], site_2[[1, 2]]])
        site_diff = sorted_sites[1] - sorted_sites[2]
        println(site_diff)

        if [site_1[3], site_2[3]] == [1, 3]
            append!(inter_dimer, [(2 * sorted_bond[2], 2 * sorted_bond[1]), (2 * sorted_bond[2], 2 * sorted_bond[1] + 1)])
        else
            append!(inter_dimer, [(2 * sorted_bond[1], 2 * sorted_bond[2]), (2 * sorted_bond[1], 2 * sorted_bond[2] + 1)])
        end
    end

    println(NLCE.nv(cluster))
    println(NLCE.edge_list(cluster))
    println(NLCE.labels(cluster))
    println(NLCE.all_coordinates(cluster))
    println(intra_dimer)
    println(inter_dimer)
    println("----------------------------")

end

end
