"""
    AbstractExpansion

Abstract base type for an NLCE expansion. Stores the subgraph relationships between
clusters and the weight matrix that `summation!` populates.

Subtypes must implement:
- `Base.getindex(e, cluster_id, order)` — weight for a cluster at a given order
- `Base.length(e)` — total number of clusters (including single-site terms)
- `add_array!(e, order, per_cluster)` — accumulate a per-cluster vector into the weights at `order`
- `order_ids(e, order)` — indices of all clusters first appearing at `order`
- `get_subclusters(e, cluster_id)` — indices of all strict subclusters of `cluster_id`
- `order_offset(e)` — `0` for site expansions, `1` for cluster expansions
"""
abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_id::Int, order::Int) = _NI("Base.getindex")
Base.length(e::AbstractExpansion) = _NI("Base.length")
add_array!(e::AbstractExpansion, order::Int, per_cluster::AbstractVector{Float64}) = _NI("add_array!")
order_ids(e::AbstractExpansion, order::Int) = _NI("order_ids")
get_subclusters(e::AbstractExpansion, cluster_id::Int) = _NI("get_subclusters")
order_offset(e::AbstractExpansion) = _NI("order_offset")

"""
    summation!(expansion, max_order)

Performs the recursive NLCE summation up till the max_order and populates the given expansion with the resultant weights
"""
summation!(e::AbstractExpansion, max_order::Int) = _summation!(e, max_order + order_offset(e))

function _summation!(e::AbstractExpansion, max_order::Int)
        stack = Vector{Tuple{Int,Float64}}()
        per_cluster = zeros(Float64, length(e))
        for order in 1:max_order
                for cluster_id in order_ids(e, order)
                        updated_lattice_constant = e[cluster_id, order]
                        push!(stack, (cluster_id, -updated_lattice_constant))
                        while !isempty(stack)
                                current_cluster_id, contribution = pop!(stack)
                                for subcluster_id in get_subclusters(e, current_cluster_id)
                                        @inbounds per_cluster[subcluster_id] += contribution
                                        push!(stack, (subcluster_id, -contribution))
                                end
                        end
                end
                add_array!(e, order, per_cluster)
                fill!(per_cluster, zero(Float64))
        end
        e
end

include("util.jl")
include("Expansion.jl")
