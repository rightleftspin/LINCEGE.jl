struct ExpansionVertex{C<:AbstractArray,B<:Bool,R<:RealSpaceVertices,L<:Integer,E<:ExpansionVertices}
        coord::C
        center::B
        real_space_vertices::R
end

center(v::ExpansionVertex) = v.center

real_space_vertices(v::ExpansionVertex) = v.real_space_vertices

function real_space_vertices(vs::AbstractVector{<:ExpansionVertex})
        foldl(union, vs)
end

function real_space_vertices(vs::AbstractSet{<:ExpansionVertex})
        foldl(union, vs)
end
