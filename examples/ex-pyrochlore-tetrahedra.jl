
using NLCE
using JLD

pyrochlore_lattice = Dict("Expansion Basis" => [[-1/2, -1/2, -1/2], [1/2, 1/2, 1/2]],
                      "Struct Per Basis" => [[[1/2, 1/2, 1/2], [1/2, -1/2, -1/2], [-1/2, 1/2, -1/2], [-1/2, -1/2, 1/2]],
                                             [[-1/2, -1/2, -1/2], [-1/2, 1/2, 1/2], [1/2, -1/2, 1/2], [1/2, 1/2, -1/2]],
                                            ],
                      "Expansion Labels" => [[1, 2, 3, 4], [1, 2, 3, 4]],
                      "Expansion Primitive Vectors" => [[2, 2, 0], [2, 0, 2], [0, 2, 2]],
                      "Expansion Neighbors" => [sqrt(3)],
                      )

# Order starts a 0 for the single site, then goes up from there
neighbors = [sqrt(2)]
max_order = 6

pyrochlore_nlce_bundle = NLCE.WeakClusterExpansionBundle(
    pyrochlore_lattice["Expansion Basis"],
    pyrochlore_lattice["Struct Per Basis"],
    pyrochlore_lattice["Expansion Labels"],
    pyrochlore_lattice["Expansion Primitive Vectors"],
    pyrochlore_lattice["Expansion Neighbors"],
    neighbors,
    max_order,
    NLCE.isomorphic_pruning,
)

pyrochlore_lattice_cluster_info = NLCE.lattice_constants!(
    pyrochlore_nlce_bundle,
    (length(unique(Iterators.flatten(pyrochlore_lattice["Expansion Labels"])))),
    single_site=true,
)

NLCE.subclusters!(pyrochlore_nlce_bundle, true)

final_weights = NLCE.final_clusters(pyrochlore_nlce_bundle, true)

num_sites = []
bond_lists = []
multiplicities = []
for (cluster, mults) in final_weights
    push!(num_sites, NLCE.nv(cluster))
    push!(bond_lists, NLCE.weighted_edge_list(cluster))
    push!(multiplicities, mults)
end

save("./outputs/Pyrochlore_Lattice_$(max_order).jld",
     "num_sites", num_sites,
     "bond_lists", bond_lists,
     "multiplicities", multiplicities
     )
