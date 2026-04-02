"""
    AbstractVertices{V}

Abstract base type for a set of vertices of element type `V`.

Subtypes must implement:
- `vertices(vs)` — return the underlying iterable of vertex indices
- `Base.collect(vs)` — return a sorted `Vector{V}`
- `Base.sort(vs)` — return a sorted copy
- `Base.intersect(vs1, vs2)` — set intersection
- `Base.setdiff(vs1, vs2)` — set difference
- `Base.union(vs1, vs2)` — set union (two-argument form)
- `Base.in(v, vs)` — membership test
- `Base.eltype(vs)` — element type `V`
"""
abstract type AbstractVertices{V} end

vertices(vs::AbstractVertices) = _NI("vertices")
Base.collect(vs::AbstractVertices) = _NI("Base.collect")
Base.sort(vs::AbstractVertices) = _NI("Base.sort")
Base.intersect(vs1::AbstractVertices, vs2::AbstractVertices) = _NI("Base.intersect")
Base.setdiff(vs1::AbstractVertices, vs2::AbstractVertices) = _NI("Base.setdiff")
Base.union(vs1::AbstractVertices, vs2::AbstractVertices) = _NI("Base.union")
Base.in(v, vs::AbstractVertices) = _NI("Base.in")
Base.eltype(vs::AbstractVertices) = _NI("Base.eltype")

Base.length(vs::AbstractVertices) = length(vertices(vs))
Base.isempty(vs::AbstractVertices) = isempty(vertices(vs))

Base.isequal(vs1::AbstractVertices, vs2::AbstractVertices) = vs1 == vs2
Base.:(==)(vs1::AbstractVertices, vs2::AbstractVertices) = (collect(vs1) == collect(vs2))
Base.hash(vs::AbstractVertices, h::UInt) = hash(collect(sort(vs)), h)
Base.contains(vs::AbstractVertices, v) = v in vs
Base.haskey(vs::AbstractVertices, v) = v in vs

function Base.union(vs::AbstractVertices, itr)
        vstemp = vs
        for x in itr
                vstemp = union(vstemp, x)
        end
        vstemp
end

Base.getindex(vec::Vector, vs::AbstractVertices) = vec[collect(vs)]
Base.getindex(mat::Matrix, vs1::AbstractVertices, vs2::AbstractVertices) = mat[collect(vs1), collect(vs2)]
Base.getindex(mat::Matrix, ::Colon, vs::AbstractVertices) = mat[:, collect(vs)]
Base.getindex(mat::Matrix, vs::AbstractVertices, ::Colon) = mat[collect(vs), :]

Base.iterate(vs::AbstractVertices) = iterate(vertices(vs))
Base.iterate(vs::AbstractVertices, state) = iterate(vertices(vs), state)

Base.show(io::IO, vs::AbstractVertices) = print(io, "Vertices: ", collect(vs))

include("TaggedVertices.jl")
