struct Vertices{V<:AbstractSet{<:Integer}}
        vertices::V
end

function Vertices()
        Vertices(BitSet())
end

function Vertices(i::Int)
        Vertices(BitSet(i))
end
vertices(vs::Vertices) = vs.vertices

Base.collect(vs::Vertices) = collect(vs.vertices)
Base.length(vs::Vertices) = length(vs.vertices)
Base.getindex(vec::Vector, vs::Vertices) = vec[collect(vs)]

Base.in(v, vs::Vertices) = v in vs.vertices
Base.intersect(vs1::Vertices, vs2::Vertices) = Vertices(intersect(vertices(vs1), vertices(vs2)))
Base.union(vs1::Vertices, vs2::Vertices) = Vertices(union(vertices(vs1), vertices(vs2)))
function Base.union(vs::Vertices, itr)
        vstemp = vs
        for x in itr
                vstemp = union(vstemp, x)
        end

        vstemp
end
