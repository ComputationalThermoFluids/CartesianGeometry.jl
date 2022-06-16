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

function _integrate!(v::ArrayAbstract{1}, coor::ArrayAbstract{1},
                     f, xyz, outer)
    (x,) = xyz

    outer = CartesianIndices(outer)

    xex = zeros(Cdouble, 4)

    for mnp in outer
        (m,) = Tuple(mnp)

        v[m] = _cell_integrate!(xex, f, SVector(x[m], x[m+1]))
        coor[m] = SVector(xex[1])
    end

    v
end

function _integrate!(v::ArrayAbstract{2}, coor::ArrayAbstract{2},
                     f, xyz, outer)
    (x, y) = xyz

    outer = CartesianIndices(outer)

    xex = zeros(Cdouble, 4)

    for mnp in outer
        (m, n) = Tuple(mnp)

        v[m, n] = _cell_integrate!(xex, f, SVector(x[m], x[m+1]),
                                           SVector(y[n], y[n+1]))
        coor[m, n] = SVector(xex[1], xex[2])
    end

    v
end

function _integrate!(v::ArrayAbstract{3}, coor::ArrayAbstract{3},
                     f, xyz, outer)
    (x, y, z) = xyz

    outer = CartesianIndices(outer)

    xex = zeros(Cdouble, 4)

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        v[m, n, p] = _cell_integrate!(xex, f, SVector(x[m], x[m+1]),
                                              SVector(y[n], y[n+1]),
                                              SVector(z[p], z[p+1]))
        coor[m, n, p] = SVector(xex[1], xex[2], xex[3])
    end

    v
end

"""

!!! note

    Even if cell is empty, return full centroid coordinates.

"""
function _cell_integrate!(xex, f, x)
    t = SVector{2}(f(i) for i in x)

    vol = x[2] - x[1]

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        return vol
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        return zero(vol)
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

    vol = (x[2] - x[1]) * (y[2] - y[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        return vol
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        return zero(vol)
    end

    x0 = Cdouble.((x[1], y[1], zero(vol)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(vol)))

    vol * getcc(f, x0, h0, xex, Cint(2); nex=Cint.((1, 0)))
end

function _cell_integrate!(xex, f, x, y, z)
    t = SArray{Tuple{2,2,2}}(f(i, j, k) for i in x, j in y, k in z)

    vol = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        return vol
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        return zero(vol)
    end

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    vol * getcc(f, x0, h0, xex, Cint(3); nex=Cint.((1, 0)))
end
