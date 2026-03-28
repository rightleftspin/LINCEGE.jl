struct Expansion <: AbstractExpansion
        index_dictionary::Dict{UInt,Int}
        subgraphs::Vector{Vector{Int}}
        weights::Matrix{Float64}
        order_ids::Dict{Int,Vector{Int}}
        order_offset::Int
end

function Expansion(clusters::AbstractClusterSet, lattice::SiteExpansionLattice, max_order::Int)
        index_dictionary = Dict{UInt,Int}()
        subgraphs = Vector{Vector{Int}}()
        weights = zeros(Float64, length(clusters), max_order)
        order_ids = Dict{Int,Vector{Int}}()

        for (i, cluster) in enumerate(sort(clusters))
                if haskey(order_ids, length(cluster))
                        push!(order_ids[length(cluster)], i)
                else
                        order_ids[length(cluster)] = [i]
                end
                index_dictionary[cluster.ghash] = i
                weights[i, length(cluster)] = cluster.lc
                temp_subgraphs = Int[]
                for subgraph_evs in get_subgraphs(cluster, lattice)
                        push!(temp_subgraphs, index_dictionary[ghash(clusters, subgraph_evs)])
                end
                push!(subgraphs, temp_subgraphs)
        end

        Expansion(index_dictionary, subgraphs, weights, order_ids, 0)
end

function Expansion(clusters::AbstractClusterSet, lattice::AbstractClusterExpansionLattice, max_order::Int)
        index_dictionary = Dict{UInt,Int}()
        subgraphs = Vector{Vector{Int}}()
        order_ids = Dict{Int,Vector{Int}}(1 => [])
        weights = zeros(Float64, length(clusters), max_order + 1)

        for (i, cluster) in enumerate(sort(clusters))
                if haskey(order_ids, length(cluster) + 1)
                        push!(order_ids[length(cluster)+1], i)
                else
                        order_ids[length(cluster)+1] = [i]
                end
                index_dictionary[cluster.ghash] = i
                weights[i, length(cluster)+1] = cluster.lc
                temp_subgraphs = Int[]

                for lv in connections(lattice)[cluster.evs]
                        gh = ghash(clusters, lv)

                        if !haskey(index_dictionary, gh)
                                ind = size(weights, 1) + 1
                                index_dictionary[gh] = ind
                                weights = vcat(weights, zeros(1, max_order + 1))
                                weights[ind, 1] = 1 / n_site_colors(lattice)
                                push!(order_ids[1], ind)
                        end

                        push!(temp_subgraphs, index_dictionary[gh])
                end

                for subgraph_evs in get_subgraphs(cluster, lattice)
                        push!(temp_subgraphs, index_dictionary[ghash(clusters, subgraph_evs)])
                end
                push!(subgraphs, temp_subgraphs)
        end

        Expansion(index_dictionary, subgraphs, weights, order_ids, 1)
end

Base.getindex(e::Expansion, cluster_id::Int, order::Int) = getindex(e.weights, cluster_id, order)
Base.length(e::Expansion) = length(e.subgraphs)
order_ids(e::Expansion, order::Int) = e.order_ids[order]
get_subclusters(e::Expansion, cluster_id::Int) = e.subgraphs[cluster_id]
add_array!(e::Expansion, order::Int, per_cluster::AbstractVector{Float64}) = @views e.weights[:, order] .+= per_cluster
order_offset(e::Expansion) = e.order_offset

function write_to_json(e::Expansion, lattice::SiteExpansionLattice, cs::AbstractClusterSet, filepath::String)
        all_coords = get_coordinates(lattice)
        all_colors = get_site_colors(lattice)
        adj = bond_matrix(lattice)

        clusters_data = []
        for cluster in cs
                n = length(cluster.vs)
                cluster_id = e.index_dictionary[cluster.ghash]

                coords = collect(eachcol(all_coords[:, cluster.vs]))
                colors = all_colors[cluster.vs]
                bonds = adj_mat_to_edge_list(adj[cluster.vs, cluster.vs])

                push!(clusters_data, Dict(
                        "cluster_id" => cluster_id,
                        "n_sites" => n,
                        "coordinates" => coords,
                        "site_colors" => colors,
                        "bonds" => bonds,
                        "weights" => vec(e.weights[cluster_id, :])
                ))
        end

        open(filepath, "w") do io
                JSON3.pretty(io, clusters_data)
        end
end
