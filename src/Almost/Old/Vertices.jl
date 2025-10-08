abstract type AbstractVertex end

dimension(v::AbstractVertex) = length(v.coord)
neighbors(v::AbstractVertex) = v.neighbors
coord(v::AbstractVertex) = v.coord

Base.hash(v::AbstractVertex, h::UInt) = hash(v.coord, h)
Base.isequal(v1::AbstractVertex, v2::AbstractVertex) = (v1.coord == v2.coord)

abstract type AbstractVertices end

vertices(vs::AbstractVertices) = vs.vertices

function Base.union(vs::AbstractVertices, itr)
        vstemp = vs
        for x in itr
                vstemp = union(vstemp, x)
        end

        vstemp
end

Base.collect(vs::AbstractVertices) = collect(vs.vertices)
Base.in(v, vs::AbstractVertices) = v in vs.vertices

Base.length(vs::AbstractVertices) = length(vs.vertices)

Base.getindex(vec::AbstractVector, vs::AbstractVertices) = vec[collect(vs)]

