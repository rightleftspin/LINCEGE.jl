abstract type AbstractCluster end

expansion_vertices(cluster::AbstractCluster) = cluster.expansion_vertices

Base.length(cluster::AbstractCluster) = length(cluster.expansion_vertices)

Base.hash(cluster::AbstractCluster, h::UInt) = hash(cluster.ghash, h)
Base.isequal(g1::AbstractCluster, g2::AbstractCluster) = (g1.ghash == g2.ghash)
