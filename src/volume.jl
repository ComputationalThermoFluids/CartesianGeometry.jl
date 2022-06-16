function integrate(f, xyz, ranges)
    (; outer) = ranges
    v = Vector{Float64}(undef, prod(length.(outer)))
    coor = Vector{SVector{length(outer),Float64}}(undef, prod(length.(outer)))
    integrate!(v, coor, f, xyz, ranges)
end

function integrate!(v, coor, f, xyz, ranges)
    (; outer, inner) = ranges

    rv = reshape(v, length.(outer)...)
    rcoor = reshape(coor, length.(outer)...)

    xyz = map(xyz, outer) do x, rng
        reshape(x, length(rng))
    end

    outer = findin.(inner, outer)

    _integrate!(rv, rcoor, f, xyz, outer)

    v, coor
end

@generated function _integrate!(v::ArrayAbstract{N}, coor::ArrayAbstract{N},
                                f, xyz, outer) where {N}
    quote
        @ntuple($N, x) = xyz

        outer = CartesianIndices(outer)

        xex = zeros(Cdouble, 4)

        for index in outer
            @ntuple($N, i) = Tuple(index)
            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))

            @nref($N, v, i) = @ncall($N, _cell_integrate!, xex, f, y)
            @nref($N, coor, i) = @ncall($N, SVector, d -> xex[d])
        end

        v
    end
end

"""

!!! note

    Even if cell is empty, return full centroid coordinates.

"""
function _cell_integrate!(xex, f, x)
    t = SVector{2}(f(i) for i in x)

    val = x[2] - x[1]

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        return zero(val)
    end

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    if isnonnegative(t[1])
        xex[1] = (ξ + x[2]) / 2
        return x[2] - ξ
    end

    if isnonnegative(t[2])
        xex[1] = (x[1] + ξ) / 2
        return ξ - x[1]
    end

    nothing
end

function _cell_integrate!(xex, f, x, y)
    t = SMatrix{2,2}(f(i, j) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    val * getcc(f, x0, h0, xex, Cint(2); nex=Cint.((1, 0)))
end

function _cell_integrate!(xex, f, x, y, z)
    t = SArray{Tuple{2,2,2}}(f(i, j, k) for i in x, j in y, k in z)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    val * getcc(f, x0, h0, xex, Cint(3); nex=Cint.((1, 0)))
end
