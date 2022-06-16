struct Dirichlet{F,T,R}
    f::F
    xyz::T
    ranges::R

    Dirichlet(f::F, xyz::T, ranges::R) where {F,T,R} =
        new{F,T,R}(f, xyz, ranges)
    Dirichlet(f::Type{F}, xyz::T, ranges::R) where {F,T,R} =
        new{Type{F},T,R}(f, xyz, ranges)
end

function apply!(v, border::Dirichlet)
    (; f, xyz, ranges) = border
    (; outer, inner) = ranges

    rv = reshape(v, length.(outer)...)

    rxyz = map(xyz, outer) do x, r
        reshape(x, length(r))
    end

    inner = CartesianGeometry.findin.(inner, outer)
    outer = eachindex.(outer)

    _apply!(rv, f, rxyz, outer, inner)

    v
end

@generated function _apply!(v::ArrayAbstract{N}, f, xyz, outer, inner) where {N}
    quote
        @ntuple($N, x) = xyz

        for index in EdgeIterator(outer, inner)
            @ntuple($N, i) = Tuple(index)
            @nref($N, v, i) = @ncall($N, f, d -> x_d[i_d])
        end

        v
    end
end
