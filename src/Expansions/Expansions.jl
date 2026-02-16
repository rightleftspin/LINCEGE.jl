module Expansions

import DataStructures.Deque

import LINCEGE:
    Vertices.ExpansionVertices,
    Lattices.AbstractLattice,
    Lattices.SiteExpansionLattice,
    Clusters.AbstractCluster,
    Clusters.AbstractClusterSet,
    Clusters.ghash

abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_id::Int, order::Int) = _NI("Base.getindex")
Base.setindex!(e::AbstractExpansion, l::Rational, cluster_id::Int, order::Int) = _NI("Base.setindex!")
order_ids(e::AbstractExpansion, order::Int) = _NI("num_clusters")

function summation!(e::AbstractExpansion, max_order::Int)
    for order in 1:max_order
        for cluster_id in order_ids(e, order)
            updated_lattice_constant = e[cluster_id]
            e[cluster_id, order] += updated_lattice_constant

            subgraphs = Deque{(Int, Int)}()
            pushfirst!(subgraphs, (cluster_id, 0))
            while !(length(subgraphs) == 0)
                current_cluster_id, current_depth = popfirst!(subgraphs)
                for subcluster_id in get_subclusters(cs, current_cluster_id)
                    e[subcluster_id, order] += ((-1)^(current_depth + 1)) * updated_lattice_constant
                    pushfirst!(subgraphs, (subcluster_id, current_depth + 1))
                end
            end
        end
    end
    e
end

include("util.jl")
include("SiteExpansions.jl")
end
