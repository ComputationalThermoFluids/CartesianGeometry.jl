function integrate(f, xyz, ranges)
    (; outer) = ranges
    v = Vector{Float64}(undef, prod(length.(outer)))
    integrate!(v, f, xyz, ranges)
end

function integrate!(v, f, xyz, ranges)
    (; outer, inner) = ranges

    w = reshape(v, length.(outer)...)

    xyz = map(xyz, outer) do x, rng
        reshape(x, length(rng))
    end

    ranges = map(ranges) do el
        findin.(inner, el)
    end

    _integrate!(w, f, xyz, ranges)

    v
end

function _integrate(f, x)
    t = f.(x) 

    all(ispositive, t) && return x[2] - x[1]
    all(isnegative, t) && return zero(eltype(x))

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    ispositive(t[1]) && return ξ - x[1]
    ispositive(t[2]) && return x[2] - ξ

    nothing
end

function _integrate!(v::ArrayAbstract{1}, f, xyz, ranges)
    (x,) = xyz

    (; outer, inner) = map(ranges) do el
        CartesianIndices(el)
    end

    for (ijk, mnp) in zip(inner, outer)
        (i,) = Tuple(ijk)
        (m,) = Tuple(mnp)

        v[m] = _integrate(f, (x[m], x[m+1]))
    end

    v
end
