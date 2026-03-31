abstract type AbstractHasher end

ghash(h::AbstractHasher, evs::ExpansionVertices) = _NI("ghash")
ghash(h::AbstractHasher, lvs::LatticeVertices) = _NI("ghash")

include("util.jl")
include("TranslationHasher.jl")
include("IsomorphicHasher.jl")
include("SymmetricHasher.jl")
