function getsurface(f, xyz, coor, ranges)
    (; outer) = ranges
    b = Vector{Float64}(undef, length(outer) * prod(length.(outer)))
    getsurface!(b, f, xyz, coor, ranges)
end

function getsurface!(b, f, xyz, coor, ranges)
    (; outer, inner) = ranges

    rb = reshape(b, length.(outer)..., :)
    rcoor = reshape(coor, length.(outer)...)

    xyz = map(xyz, outer) do x, rng
        reshape(x, length(rng))
    end

    outer = findin.(inner, outer)

    _getsurface!(rb, f, xyz, rcoor, outer)

    b
end

function _getsurface!(b::ArrayAbstract{2}, f, xyz, coor::ArrayAbstract{1}, outer)
    (x,) = xyz

    outer = CartesianIndices(outer)

    xex = zeros(Cdouble, 4)

    # x-normal faces
    for mnp in outer
        (m,) = Tuple(mnp)

        b[m, 1] = _getpoint!(xex, f, coor[m][1])
    end

    b
end

function _getsurface!(b::ArrayAbstract{3}, f, xyz, coor::ArrayAbstract{2}, outer)
    (x, y) = xyz

    outer = CartesianIndices(outer)

    xex = zeros(Cdouble, 4)

    # x-normal faces
    for mnp in outer
        (m, n) = Tuple(mnp)

        b[m, n, 1] = _getlength!(xex, f, coor[m, n][1], SVector(y[n], y[n+1]))
    end

    # y-normal faces
    for mnp in outer
        (m, n) = Tuple(mnp)

        b[m, n, 2] = _getlength!(xex, f, SVector(x[m], x[m+1]), coor[m, n][2])
    end

    b
end

function _getsurface!(b::ArrayAbstract{4}, f, xyz, coor::ArrayAbstract{3}, outer)
    (x, y, z) = xyz

    outer = CartesianIndices(outer)

    xex = zeros(Cdouble, 4)

    # x-normal faces
    for mnp in outer
        (m, n, p) = Tuple(mnp)

        b[m, n, p, 1] = _getarea!(xex, f, coor[m, n, p][1], SVector(y[n], y[n+1]),
                                                            SVector(z[p], z[p+1]))
    end

    # y-normal faces
    for mnp in outer
        (m, n, p) = Tuple(mnp)

        b[m, n, p, 2] = _getarea!(xex, f, SVector(x[m], x[m+1]), coor[m, n, p][2],
                                          SVector(z[p], z[p+1]))
    end

    # z-normal faces
    for mnp in outer
        (m, n, p) = Tuple(mnp)

        b[m, n, p, 3] = _getarea!(xex, f, SVector(x[m], x[m+1]),
                                          SVector(y[n], y[n+1]), coor[m, n, p][3])
    end

    b
end

