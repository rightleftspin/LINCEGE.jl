function get_subgraphs(c::AbstractCluster, lattice::SiteExpansionLattice)
    max_depth = length(c) - 1
    roots = [ExpansionVertices(center) for center in c.evs]
    println(roots)
    visited = Set(roots)
    vlock = ReentrantLock()

    function try_mark(cluster::ExpansionVertices)
        lock(vlock)
        already = cluster in visited
        if !already
            push!(visited, cluster)
        end
        unlock(vlock)

        if length(cluster) == max_depth
            return false
        end
        !already
    end

    function dfs(cluster::ExpansionVertices)
        if !try_mark(cluster)
            return
        end
        for ev in neighbors(lattice, cluster.evs)
            if ev in c.evs
                dfs(union(cluster.evs, ExpansionVertices(ev)))
            end
        end
    end

    for c in roots
        dfs(c)
    end

    visited
end
