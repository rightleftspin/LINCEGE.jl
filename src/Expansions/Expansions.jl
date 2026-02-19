module Expansions

using DataStructures

import LINCEGE:
    Vertices.ExpansionVertices,
    Lattices.AbstractLattice,
    Lattices.SiteExpansionLattice,
    Lattices.neighbors,
    Clusters.AbstractCluster,
    Clusters.AbstractClusterSet,
    Clusters.ghash

abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_id::Int, order::Int) = _NI("Base.getindex")
Base.getindex(e::AbstractExpansion, cluster_ids::Vector{Int}, order::Int) = _NI("Base.getindex")
Base.setindex!(e::AbstractExpansion, l::Rational, cluster_id::Int, order::Int) = _NI("Base.setindex!")
Base.setindex!(e::AbstractExpansion, l::Rational, cluster_ids::Vector{Int}, order::Int) = _NI("Base.setindex!")
Base.setindex!(e::AbstractExpansion, ls, cluster_ids::Vector{Int}, order::Int) = _NI("Base.setindex!")
order_ids(e::AbstractExpansion, order::Int) = _NI("num_clusters")
get_subclusters(e::AbstractExpansion, cluster_id::Int) = _NI("get_subclusters")

function summation!(e::AbstractExpansion, max_order::Int)
    subgraphs = Vector{Tuple{Int,Int}}()
    for order in 1:max_order
        for cluster_id in order_ids(e, order)
            push!(subgraphs, (cluster_id, 0))
            while !(length(subgraphs) == 0)
                current_cluster_id, current_depth = pop!(subgraphs)
                ids = get_subclusters(e, current_cluster_id)
                depths = repeat([current_depth + 1], length(ids))
                e[ids, order] += ((-1) .^ depths) .* e[cluster_id, order]
                append!(subgraphs, zip(ids, depths))
            end
        end
    end
    e
end

include("util.jl")
include("SiteExpansions.jl")
end
