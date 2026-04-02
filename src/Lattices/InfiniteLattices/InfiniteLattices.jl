"""
    AbstractInfiniteLattice <: AbstractLattice

Abstract base type for an infinite lattice used to enumerate clusters.

Subtypes must implement all methods of `AbstractLattice`, plus:
- `n_unique_sites(lattice)` — number of translationally-inequivalent sites in the lattice
"""
abstract type AbstractInfiniteLattice <: AbstractLattice end

"""
    AbstractClusterExpansionLattice <: AbstractInfiniteLattice

Abstract base type for lattices where each expansion vertex represents a cluster of
physical sites (a unit cell), rather than a single site.

Subtypes must implement all methods of `AbstractInfiniteLattice`, plus:
- `n_site_colors(lattice)` — number of distinct site colors in the lattice 
- `connections(lattice)` — mapping from expansion vertices to their constituent lattice vertices
"""
abstract type AbstractClusterExpansionLattice <: AbstractInfiniteLattice end

n_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")
n_site_colors(lattice::AbstractClusterExpansionLattice) = _NI("n_site_colors")
connections(lattice::AbstractClusterExpansionLattice) = _NI("connections")

include("util.jl")
include("SiteExpansionLattices.jl")
include("StrongClusterExpansionLattices.jl")
include("WeakClusterExpansionLattices.jl")
