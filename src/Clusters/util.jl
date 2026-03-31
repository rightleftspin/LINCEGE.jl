function add_cluster!(cs::AbstractClusterSet{C,H}, c::AbstractCluster) where {C,H}
        new_cluster = C(c.vs, c.lc, ghash(cs, c))
        if new_cluster in cs
                old_cluster = pop!(cs, new_cluster)
                push!(cs, C(c.vs, c.lc + old_cluster.lc, old_cluster.ghash))
        else
                push!(cs, new_cluster)
        end
end

"""
    clusters_from_clusters!(new_clusters, old_clusters)

# TODO: describe what clusters_from_clusters! does
"""
function clusters_from_clusters!(new_clusters::AbstractClusterSet{C,H}, old_clusters::AbstractClusterSet) where {C<:AbstractCluster,H}
        for old_cluster in old_clusters
                add_cluster!(new_clusters, old_cluster)
        end
        new_clusters
end

"""
    clusters_from_lattice!(clusters, lattice; spawn_depth=3)

# TODO: describe what clusters_from_lattice! does
"""
function clusters_from_lattice!(clusters::AbstractClusterSet{C,H}, lattice::AbstractInfiniteLattice; spawn_depth::Int=3) where {C<:AbstractCluster,H}
        max_depth = max_order(lattice)
        ctrs = centers(lattice)
        roots = [C(typeof(ctrs)(center), clusters, lattice) for center in ctrs]
        vlock = ReentrantLock()

        function try_mark(cluster::AbstractCluster)
                already = lock(vlock) do
                        in_cs = cluster in clusters
                        if !in_cs
                                push!(clusters, cluster)
                        end
                        in_cs
                end

                if length(cluster) == max_depth
                        return false
                end
                !already
        end

        function dfs(cluster::AbstractCluster, depth)
                if !try_mark(cluster)
                        return
                end

                if depth < spawn_depth
                        for v in neighbors(lattice, cluster.vs)
                                dfs(C(union(cluster.vs, typeof(ctrs)(v)), clusters, lattice), depth + 1)
                        end
                else
                        tasks = Task[]
                        first = true
                        for v in neighbors(lattice, cluster.vs)
                                neighbor_cluster = C(union(cluster.vs, typeof(ctrs)(v)), clusters, lattice)

                                if first
                                        dfs(neighbor_cluster, depth + 1)
                                        first = false
                                else
                                        push!(tasks, @spawn dfs(neighbor_cluster, depth + 1))
                                end
                        end
                        for t in tasks
                                fetch(t)
                        end
                end
        end

        for c in roots
                dfs(c, 0)
        end

        clusters
end
