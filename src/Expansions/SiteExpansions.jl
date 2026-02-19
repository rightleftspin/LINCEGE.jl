struct SiteExpansion <: AbstractExpansion
    index_dictionary::AbstractDict{UInt,Int}
    subgraphs::AbstractVector{<:AbstractVector{Int}}
    weights::AbstractMatrix{Rational}
    order_ids::AbstractDict{Int,<:AbstractVector{Int}}
end

function SiteExpansion(clusters::AbstractClusterSet, lattice::SiteExpansionLattice, max_order::Int)
    index_dictionary = Dict{UInt,Int}()
    subgraphs = Vector{Vector{Int}}()
    weights = zeros(Rational, length(clusters), max_order)
    order_ids = Dict{Int,Vector{Int}}()
    for (i, cluster) in enumerate(sort(clusters))
        if haskey(order_ids, length(cluster))
            push!(order_ids[length(cluster)], i)
        else
            order_ids[length(cluster)] = [i]
        end
        index_dictionary[cluster.ghash] = i
        weights[i, length(cluster)] = cluster.lc
        temp_subgraphs::Vector{Int} = []
        for subgraph_evs in get_subgraphs(cluster, lattice)
            push!(temp_subgraphs, index_dictionary[ghash(clusters, subgraph_evs)])
        end
        push!(subgraphs, temp_subgraphs)
    end

    SiteExpansion(
        index_dictionary,
        subgraphs,
        weights,
        order_ids
    )
end

Base.getindex(e::SiteExpansion, cluster_id::Int, order::Int) = getindex(e.weights, cluster_id, order)
Base.getindex(e::SiteExpansion, cluster_ids::Vector{Int}, order::Int) = getindex(e.weights, cluster_ids, order)
Base.setindex!(e::SiteExpansion, l::Rational, cluster_id::Int, order::Int) = setindex!(e.weights, l, cluster_id, order)
Base.setindex!(e::SiteExpansion, l::Rational, cluster_ids::Vector{Int}, order::Int) = setindex!(e.weights, l, cluster_ids, order)
Base.setindex!(e::SiteExpansion, ls, cluster_ids::Vector{Int}, order::Int) = setindex!(e.weights, ls, cluster_ids, order)
order_ids(e::SiteExpansion, order::Int) = e.order_ids[order]
get_subclusters(e::SiteExpansion, cluster_id::Int) = e.subgraphs[cluster_id]
