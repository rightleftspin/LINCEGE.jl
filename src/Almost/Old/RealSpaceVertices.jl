struct RealSpaceVertices{V<:AbstractSet{<:Integer}}
        vertices::V
end

Base.union(vs1::RealSpaceVertices, vs2::RealSpaceVertices) = RealSpaceVertices(union(vertices(vs1), vertices(vs2)))

Base.intersect(vs1::RealSpaceVertices, vs2::RealSpaceVertices) = RealSpaceVertices(intersect(vertices(vs1), vertices(vs2)))
