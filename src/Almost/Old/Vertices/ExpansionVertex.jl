struct ExpansionVertex{C<:AbstractArray,B<:Bool,R<:RealSpaceVertices,L<:Integer,E<:ExpansionVertices}
        coord::C
        center::B
end

center(v::ExpansionVertex) = v.center
