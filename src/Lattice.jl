
abstract type AbstractLattice end

mutable struct Lattice{C,V,N,W,F} <: AbstractLattice

        "Coordinates for the lattice, in sorted order"
        coordinates::AbstractMatrix{C}
        "Bonds between sites in the super lattice. Here, the index integers are instead in reference to bundles in the coordinate bundles"
        adj_list::AbstractVector{<:AbstractVector{<:Integer}}
        "Connections between the adjacency list and adjacency matrix"
        connections::AbstractVector{N}
        "Bonds between sites where the first index is a choice of representation of the cluster, second index is site 1 and third index is site 2.
        The first two representation choices are reserved for isomorphic and translational hashing respectively. The third rep. choice is reserved
        for connecting bonds to sites in the super lattice.
        Since there are no self loops, the diagonal of the first adjacency matrix contains the labels of each site, the diagonals of the
        second adjacency matrix contains the labels for each site that is used for translational invariance"
        adj_matrices::AbstractArray{<:Integer,3}

        "Start points on the lattice super vertices for the NLCE expansion"
        start::AbstractVector{V}
        "Maximum order for the NLCE expansion"
        max_order::Integer
        "Hashing function for lattice embedding constant count"
        hashing_fxn::F
        #function Lattice(
        #        lattice::Cluster,
        #        coordinates::AbstractMatrix{<:Real},
        #        start::AbstractVector{<:Integer},
        #        max_order::Integer,
        #        hashing_fxn,
        #)

        #        @assert nv(lattice) == size(coordinates, 1) "Number of coordinates needs to be the number of sites in the lattice"
        #        @assert issorted(eachrow(coordinates)) "Coordinates must be sorted"
        #        @assert maximum(start) <= nsv(lattice) "Some start points are not within the lattice"

        #        bundle = new()
        #        bundle.lattice = lattice
        #        bundle.coordinates = coordinates
        #        bundle.start = start
        #        bundle.max_order = max_order
        #        bundle.hashing_fxn = hashing_fxn
        #        bundle
        #end
end
