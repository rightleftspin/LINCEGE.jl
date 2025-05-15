
abstract type AbstractBundle end

# TODO: Documentation
"""
"""
mutable struct SiteExpansionBundle <: AbstractBundle

    "Underlying lattice for this NLCE expansion"
    lattice::Cluster
    "Coordinates for the lattice, in sorted order"
    coordinates::AbstractMatrix{<:Real}
    "Sublattice Coordinates for the lattice, in sorted order"
    sublattice_coordinates::AbstractMatrix{<:Integer}
    "Start points on the lattice super vertices for the NLCE expansion"
    start::AbstractVector{<:Integer}
    "Maximum order for the NLCE expansion"
    max_order::Integer
    "Hashing function for lattice embedding constant count"
    hashing_fxn
    "Cluster info for all clusters invariant under the hashing function of this bundle"
    cluster_info

    function SiteExpansionBundle(
        lattice::Cluster,
        coordinates::AbstractMatrix{<:Real},
        sublattice_coordinates::AbstractMatrix{<:Integer},
        start::AbstractVector{<:Integer},
        max_order::Integer,
        hashing_fxn
    )

        @assert nv(lattice) == nsv(lattice) "Number of Vertices and Super Vertices need to be equivalent for site expansion"
        @assert nv(lattice) == size(coordinates, 1) "Number of coordinates needs to be the number of sites in the lattice"
        @assert nv(lattice) == size(sublattice_coordinates, 1) "Number of sublattice coordinates needs to be the number of sites in the lattice"
        @assert size(coordinates, 2) == (size(sublattice_coordinates, 2) - 1) "Sublattice coordinate or real space coordiante dimension is incorrect"
        @assert issorted(eachrow(coordinates)) "Coordinates must be sorted"
        @assert maximum(start) <= nsv(lattice) "Some start points are not within the lattice"

        bundle = new()
        bundle.lattice = lattice
        bundle.coordinates = coordinates
        bundle.sublattice_coordinates = sublattice_coordinates
        bundle.start = start
        bundle.max_order = max_order
        bundle.hashing_fxn = hashing_fxn
        bundle
    end
end

"""
Generates SiteExpansionBundle with lattice of size 2 * max_order + 1
"""
function SiteExpansionBundle(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighbors::AbstractVector{<:Real},
    max_order::Integer,
    hashing_fxn;
    basis_labels::AbstractVector{<:Integer} = repeat([1], length(basis)),
)
    # TODO: Move to util functions
    primitive_lattice, unrotated_lattice = generate_primitive_lattice(primitive_vectors, max_order)

    unsorted_coords = add_basis_coords(basis, primitive_lattice)
    unsorted_sublattice_coords = add_basis_sublattice(basis, unrotated_lattice)

    sort_perm_coords = sortperm(eachrow(unsorted_coords))

    coords = unsorted_coords[sort_perm_coords, :]
    sublattice_coords = unsorted_sublattice_coords[sort_perm_coords, :]
    labels = repeat(basis_labels, fld(size(unsorted_coords, 1), length(basis)))[sort_perm_coords]
    translation_labels = repeat(1:length(basis), fld(size(unsorted_coords, 1), length(basis)))[sort_perm_coords]

    adj_list = adj_list_from_coords(coords, neighbors)
    adj_matrices = adj_matrices_strong(coords,
                                       neighbors,
                                       labels,
                                       translation_labels,
                                       1:size(unsorted_coords, 1)
                                       )

    start = findfirst.(isapprox.(basis), (eachrow(coords), ))

    lattice = Cluster(
        adj_list,
        [[i] for i in 1:size(coords, 1)],
        adj_matrices,
        length(unique(basis_labels)) != 1,
        length(neighbors) != 1
    )

    SiteExpansionBundle(
        lattice,
        coords,
        sublattice_coords,
        start,
        max_order,
        hashing_fxn
    )

end

begin # Access Methods
    lattice(bundle::AbstractBundle) = bundle.lattice
    start(bundle::AbstractBundle) = bundle.start
    max_order(bundle::AbstractBundle) = bundle.max_order
    dimensions(bundle::AbstractBundle) = size(bundle.coordinates, 2)
    hashing_fxn(bundle::AbstractBundle) = bundle.hashing_fxn
    cluster_info(bundle::AbstractBundle) = bundle.cluster_info
end

begin # Writer Methods
    function set_cluster_info!(bundle::AbstractBundle, cluster_info)
        bundle.cluster_info = cluster_info

        cluster_info
    end
end

begin # Individual Cluster Methods
    find_cluster(bundle::AbstractBundle, cluster::Cluster) = findfirst(cluster, bundle.clusters)
    hash_cluster(bundle::AbstractBundle, cluster::Cluster) = bundle.hashing_fxn(cluster)

    lattice_constant(bundle::SiteExpansionBundle, cluster::Cluster) = bundle.lattice_constants[find_cluster(bundle, cluster)]
    # would be nice to be able to plot a cluster and plot a lattice here
end

begin # NLCE methods
    function lattice_constants!(bundle::SiteExpansionBundle)

        t_i_clusters, super_vertices = translationally_invariant_clusters(lattice(bundle), start(bundle), max_order(bundle))

        cluster_info = lattice_constants(hashing_fxn(bundle), length(start(bundle)), t_i_clusters, super_vertices)

        set_cluster_info!(bundle, cluster_info)
    end

    function subclusters!(bundle::SiteExpansionBundle; hash_fxn = hashing_fxn(bundle))

        for (hash, (cluster, mult, perm, svs, subclusters)) in cluster_info(bundle)
            cluster_info(bundle)[hash] = (cluster, mult, perm, svs, lattice_constants(hash_fxn, 1, find_subclusters(cluster, false)...))
        end

        cluster_info(bundle)
    end

    function final_clusters!(bundle::SiteExpansionBundle)
        # Initialize an empty output dictionary
        output_dict = Dict{Cluster, Vector{<:Real}}()

        # Return the final sum for all clusters
        for order = 1:max_order(bundle)
            for (hash, mult) in nlce_summation(cluster_info(bundle), order)
                output_dict[cluster_info(bundle)[hash][1]] =
                    append!(get(output_dict, cluster_info(bundle)[hash][1], Vector{Real}()), mult)
            end
        end

        output_dict
    end

end
