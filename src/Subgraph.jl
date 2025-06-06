struct Subgraph{V<:Unsigned,H<:Unsigned}
        vertices::AbstractSet{T} # needs to be a sorted container since it is in order
        neighbors::AbstractSet{T}
        translational_hash::H
end

function Subgraph(vertex::T, lattice::AbstractLattice) where {T<:Unsigned}
        Subgraph{T}(BitSet(vertex), neighbors(lattice, vertex), translational_hash(lattice, vertex))
end

function Subgraph(vertices::AbstractSet{T}, neighbors::AbstractSet{T}, lattice::AbstractLattice) where {T<:Unsigned}
        Subgraph{T}(vertices, setdiff(neighbors, vertices), translational_hash(lattice, vertices))
end

function neighbor_subgraphs(subgraph::Subgraph, lattice::AbstractLattice)
        neighbor_subgraph.(subgraph, lattice, subgraph.neighbors)
end

function neighbor_subgraph(subgraph::Subgraph, lattice::AbstractLattice, n::Unsigned)
        Subgraph(union(subgraph.vertices, n), union(subgraph.neighbors, neighbors(lattice, n)), lattice)
end

nodes(subgraph::Subgraph) = subgraph.vertices

Base.length(subgraph::Subgraph) = length(subgraph.vertices)

Base.hash(subgraph::Subgraph, h::UInt) = hash(subgraph.translational_hash, h)
Base.isequal(g1::Subgraph, g2::Subgraph) = (g1.translational_hash == g2.translational_hash)

