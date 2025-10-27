struct Subgraphs{H<:VertexHash,C<:Subgraph} <: AbstractClusters{H,C}
    clusters::Dict{VertexHash,C}
end

function Subgraphs()
    Subgraphs{VertexHash,Subgraph}(Dict{VertexHash,Subgraph}())
end

function Subgraphs(vs::AbstractVertices, cluster::AbstractCluster, lattice::AbstractLattice)
    subgraphs = Subgraphs()
    for v in vs
        subgraph = Subgraph(typeof(vs)(v), cluster, lattice)
        subgraphs[ghash(subgraph)] = subgraph
    end

    subgraphs
end

function Subgraphs(cluster::AbstractCluster, lattice::AbstractSiteExpansionLattice)
    if length(cluster) == 1
        return Subgraphs()
    end

    max_order = length(cluster) - 1
    roots = Subgraphs(vertices(cluster), cluster, lattice)
    visited = Subgraphs()

    function try_mark(subgraph::Subgraph)
        already = subgraph in visited
        if !already
            visited[ghash(subgraph)] = subgraph
        end

        if length(subgraph) == max_order
            return false
        end
        !already
    end

    function dfs(subgraph::Subgraph)
        if !try_mark(subgraph)
            return
        end

        nbrs = neighbor_subgraphs(subgraph, lattice)
        for ssg in nbrs
            dfs(ssg)
        end
    end

    for (ghash, sg) in roots
        dfs(sg)
    end

    visited
end

_clusters(cs::Subgraphs) = cs.clusters
Base.length(cs::Subgraphs) = length(_clusters(cs))
Base.iterate(cs::Subgraphs) = iterate(_clusters(cs))
Base.iterate(cs::Subgraphs, state) = iterate(_clusters(cs), state)
Base.getindex(cs::Subgraphs, ghash::VertexHash) = getindex(cs.clusters, ghash)
Base.haskey(cs::Subgraphs, ghash::VertexHash) = haskey(cs.clusters, ghash)
Base.setindex!(cs::Subgraphs, cluster::Subgraph, ghash::VertexHash) = setindex!(cs.clusters, cluster, ghash)
