abstract type AbstractInfiniteLattice <: AbstractLattice end
abstract type AbstractClusterExpansionLattice <: AbstractInfiniteLattice end

n_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")

n_site_colors(lattice::AbstractClusterExpansionLattice) = _NI("n_site_colors")
connections(lattice::AbstractClusterExpansionLattice) = _NI("connections")

include("util.jl")
include("SiteExpansionLattices.jl")
include("StrongClusterExpansionLattices.jl")
