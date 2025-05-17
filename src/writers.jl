"""
These are writers that write a set of clusters to a file for future use. There are many styles
of writing to disk, but the standard will be JSON files
"""

function write_to_file_coordinates(
    nlce_output::AbstractDict{Cluster,Vector{<:Real}},
    cluster_hashes::AbstractDict{Cluster,Integer},
    cluster_perms::AbstractDict{Cluster,Vector{<:Integer}},
    filename::AbstractString,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        perm = cluster_perms[cluster]
        write(nlce_file, "$(nv(cluster)):")
        write(nlce_file, " $(cluster_hashes[cluster]):")
        for edge in edge_list(cluster)
            write(
                nlce_file,
                " $(join((findall(x -> x == edge[1], perm)[1], findall(x -> x == edge[2], perm)[1], edge[3]), ' '))",
            )
        end
        write(nlce_file, ":")
        for coord in all_coordinates(cluster)[perm]
            write(nlce_file, " ($(join(coord, ',')))")
        end
        write(nlce_file, ": $(join(orbits(cluster)[perm], ' '))")
        write(nlce_file, ": $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file(
    nlce_output::AbstractDict{<:Cluster,Vector{<:Real}},
    filename::AbstractString,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        write(nlce_file, "$(nv(cluster)):")
        for edge in edge_list(cluster)
            write(nlce_file, " $(join(edge, ' '))")
        end
        write(nlce_file, " : $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file_colors(
    nlce_output::AbstractDict{Cluster,Vector{<:Real}},
    filename::AbstractString,
    bond_info,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        write(nlce_file, "$(nv(cluster)):")
        for edge in get(bond_info, cluster, 0)[1]
            write(nlce_file, " $(join(edge, ' '))")
        end
        write(nlce_file, ": $(join(labels(cluster), ' ')):")
        for edge in get(bond_info, cluster, 0)[2]
            write(nlce_file, " $(join(edge, ' '))")
        end
        #for edge in edge_list(cluster)
        #    write(nlce_file, " $(join(edge, ' '))")
        #end
        #write(nlce_file, ": $(join(labels(cluster), ' ')):")
        #for coord in all_coordinates(cluster)
        #    write(nlce_file, " ($(join(coord, ',')))")
        #end
        #write(nlce_file, ": $(join(orbits(cluster), ' '))")
        write(nlce_file, ": $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file_fortran(
    nlce_output::AbstractDict{<:Cluster,Vector{<:Real}},
    filename::AbstractString,
    max_order::Integer,
)

    nlce_files = [open(filename * "_$(i).txt", "w") for i = 1:max_order]
    sorted_clusters = sort(collect(keys(nlce_output)), by = nv)

    for cluster in sorted_clusters
        edges = edge_list(cluster)
        write(nlce_files[nv(cluster)], "$(length(edges))\n")
        for edge in edges
            write(nlce_files[nv(cluster)], "$(join(edge, '\t'))\n")
        end
        write(nlce_files[nv(cluster)], "\n")
    end

    for cluster in sorted_clusters
        write(nlce_files[nv(cluster)], "$(join(nlce_output[cluster], ' '))\n")
    end

    close.(nlce_files)
end
