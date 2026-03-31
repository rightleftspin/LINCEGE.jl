struct ClusterSet{C<:AbstractCluster,H<:AbstractHasher} <: AbstractClusterSet{C,H}
        clusters::Set{C}
        hasher::H
end

"""
    TranslationClusterSet(lattice)

# TODO: describe what a TranslationClusterSet represents
"""
function TranslationClusterSet(lattice::AbstractInfiniteLattice)
        ClusterSet{Cluster,TranslationHasher}(
                Set{Cluster}(),
                TranslationHasher(lattice)
        )
end

"""
    IsomorphicClusterSet(lattice)

# TODO: describe what an IsomorphicClusterSet represents
"""
function IsomorphicClusterSet(lattice::AbstractLattice)
        ClusterSet{Cluster,IsomorphicHasher}(
                Set{Cluster}(),
                IsomorphicHasher(lattice)
        )
end

function SymmetricClusterSet(lattice::AbstractLattice, symmetries::Vector{Matrix{Float64}})
        ClusterSet{Cluster,SymmetricHasher}(
                Set{Cluster}(),
                SymmetricHasher(lattice, symmetries)
        )
end

Base.length(cs::ClusterSet) = length(cs.clusters)
Base.in(c::C, cs::ClusterSet{C,H}) where {C<:AbstractCluster,H} = c in cs.clusters
Base.iterate(cs::ClusterSet) = iterate(cs.clusters)
Base.iterate(cs::ClusterSet, state) = iterate(cs.clusters, state)
Base.push!(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = push!(cs.clusters, c)
Base.pop!(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = pop!(cs.clusters, c)
Base.sort(cs::ClusterSet) = sort(collect(cs.clusters), by=length)

ghash(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = ghash(cs.hasher, c.vs)
ghash(cs::ClusterSet, vs::AbstractVertices) = ghash(cs.hasher, vs)
