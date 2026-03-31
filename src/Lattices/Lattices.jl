abstract type AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers")
max_order(lattice::AbstractLattice) = _NI("max_order")
neighbors(lattice::AbstractLattice, vs::AbstractVertices) = _NI("neighbors")
get_coordinates(lattice::AbstractLattice) = _NI("get_coordinates")
get_labels(lattice::AbstractLattice) = _NI("get_labels")
get_site_colors(lattice::AbstractLattice) = _NI("get_site_colors")
bond_matrix(lattice::AbstractLattice) = _NI("bond_matrix")

include("Connections.jl")
include("InfiniteLattices/InfiniteLattices.jl")
