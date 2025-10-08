struct ExpansionVertices{V<:AbstractSet{<:Integer}}
        vertices::V
end

Base.union(vs1::ExpansionVertices, vs2::ExpansionVertices) = ExpansionVertices(union(vertices(vs1), vertices(vs2)))

Base.intersect(vs1::ExpansionVertices, vs2::ExpansionVertices) = ExpansionVertices(intersect(vertices(vs1), vertices(vs2)))
