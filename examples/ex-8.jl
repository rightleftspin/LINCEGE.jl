module ex8

using NLCE
using NautyGraphs

basis = [[0, 0], [1, 0], [0, 1], [1, 1]]

colors = [1, 2, 2, 1]

primitive_vec = [[2, 0], [0, 2]]

neighborhood = [1]

max_order = 8

nlce_clusters = NLCE.site_color_NLCE(basis, colors, primitive_vec, neighborhood, max_order)

filepath = "examples/outputs/ex-8/ssl_dimer"
mkpath(filepath)
filename = filepath * "/ssl_dimer_$(max_order)"



# Rewriting clusters to be useful
bond_info = Dict{NLCE.AbstractNLCECluster, Vector}()
for (cluster, multiplicities) in nlce_clusters
    intra_dimer = []
    for site in 0:(NLCE.nv(cluster) - 1)
        push!(intra_dimer, (2 * site, (2 * site + 1)))
    end

    inter_dimer = []
    #
    for bond in NLCE.edge_list(cluster)
        # Find the coordinates for each site in the bond
        site_1 = NLCE.all_coordinates(cluster)[bond[1]]
        site_2 = NLCE.all_coordinates(cluster)[bond[2]]
        site_1 = site_1[1] .* primitive_vec[1] +
            site_1[2] .* primitive_vec[2] + basis[site_1[3]]
        site_2 = site_2[1] .* primitive_vec[1] +
            site_2[2] .* primitive_vec[2] + basis[site_2[3]]

        # Sort them so that you are going from the left and bottom site to the right and top site
        site_perm = sortperm([site_1, site_2])
        # Figure out if they are horizontal (1) or vertical (2), choice here is arbitrary
        site_labels = NLCE.labels(cluster)[[bond[1], bond[2]]][site_perm]
        # reorder the bond and subtract one for vertex counting
        new_bond = bond[[1, 2]][site_perm] .- 1
        # Take the sites and sort them in bottom left to top right order
        sites = [site_1, site_2][site_perm]
        # find the direction vector from the bottom left site to the top right site
        site_diff = sites[2] - sites[1]

        if site_diff == [1, 0]
        # To the right
            if site_labels[1] == 1
            # Horizontal dimer
                append!(inter_dimer,
                        [(2 * new_bond[1] + 1, 2 * new_bond[2]),
                        (2 * new_bond[1] + 1, 2 * new_bond[2] + 1)])
            else
            # Vertical dimer
                append!(inter_dimer,
                        [(2 * new_bond[2], 2 * new_bond[1]),
                        (2 * new_bond[2], 2 * new_bond[1] + 1)])
            end
        else
        # Upwards
        if site_diff != [0, 1]
            println("error")
        end
            if site_labels[1] == 1
            # Horizontal dimer
                append!(inter_dimer,
                        [(2 * new_bond[2], 2 * new_bond[1]),
                        (2 * new_bond[2], 2 * new_bond[1] + 1)])
            else
            # Vertical dimer
                append!(inter_dimer,
                        [(2 * new_bond[1] + 1, 2 * new_bond[2]),
                        (2 * new_bond[1] + 1, 2 * new_bond[2] + 1)])
            end

        end
    end
   # new_cluster = NautyGraph(3 * NLCE.nv(cluster))

   # final_vertex_counter = 2 * NLCE.nv(cluster)
   #     for edge in intra_dimer
   #         NautyGraphs.add_edge!(new_cluster, edge[1] + 1, final_vertex_counter + 1)
   #         NautyGraphs.add_edge!(new_cluster, edge[2] + 1, final_vertex_counter + 1)
   #         final_vertex_counter += 1
   #     end

   #     for edge in inter_dimer
   #         NautyGraphs.add_edge!(new_cluster, (edge .+ 1)...)
   #     end
   #
    #orbit = orbits(new_cluster)[1:(2 * NLCE.nv(cluster))]


    bond_info[cluster] = [intra_dimer, inter_dimer]

end

NLCE.write_to_file_colors(nlce_clusters, filename, bond_info)

end
