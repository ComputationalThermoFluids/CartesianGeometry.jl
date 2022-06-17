function getvolume(f, xyz, coor, ranges)
    (; outer) = ranges
    w = Vector{Float64}(undef, length(outer) * prod(length.(outer)))
    getvolume!(w, f, xyz, coor, ranges)
end

function getvolume!(w, f, xyz, coor, ranges)
    (; outer, stagger, center) = ranges

    rw = reshape(w, length.(outer)..., :)
    coor = reshape(coor, length.(outer)...)

    xyz = map(xyz, outer) do x, rng
        reshape(x, length(rng))
    end

    stagger = findin.(stagger, outer)
    center = findin.(center, outer)

    _getvolume!(rw, f, xyz, coor, stagger, center)

    w
end

function _getvolume!(w::ArrayAbstract{2}, f, xyz, coor, stagger, center)
    (x,) = xyz

    xex = zeros(Cdouble, 4)

    # x-normal faces
    outer = CartesianIndices((stagger[1],))

    for mnp in outer
        (m,) = Tuple(mnp)

        w[m, 1] = _getlength!(xex, f, SVector(coor[m-1][1], coor[m][1]))
    end

    w
end

function _getvolume!(w::ArrayAbstract{3}, f, xyz, coor, stagger, center)
    (x, y) = xyz

    xex = zeros(Cdouble, 4)

    # x-normal faces
    outer = CartesianIndices((stagger[1], center[2]))

    for mnp in outer
        (m, n) = Tuple(mnp)

        w[m, n, 1] = _getarea!(xex, f, SVector(coor[m-1, n][1], coor[m, n][1]),
                                       SVector(y[n], y[n+1]))
    end

    # y-normal faces
    outer = CartesianIndices((center[1], stagger[2]))

    for mnp in outer
        (m, n) = Tuple(mnp)

        w[m, n, 2] = _getarea!(xex, f, SVector(x[m], x[m+1]),
                                       SVector(coor[m, n-1][2], coor[m, n][2]))
    end

    w
end

function _getvolume!(w::ArrayAbstract{4}, f, xyz, coor, stagger, center)
    (x, y, z) = xyz

    xex = zeros(Cdouble, 4)

    # x-normal faces
    outer = CartesianIndices((stagger[1], center[2], center[3]))

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        w[m, n, p, 1] = _getvolume!(xex, f, SVector(coor[m-1, n, p][1], coor[m, n, p][1]),
                                            SVector(y[n], y[n+1]),
                                            SVector(z[p], z[p+1]))
    end

    # y-normal faces
    outer = CartesianIndices((center[1], stagger[2], center[3]))

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        w[m, n, p, 2] = _getvolume!(xex, f, SVector(x[m], x[m+1]),
                                            SVector(coor[m, n-1, p][2], coor[m, n, p][2]),
                                            SVector(z[p], z[p+1]))
    end

    # z-normal faces
    outer = CartesianIndices((center[1], center[2], stagger[3]))

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        w[m, n, p, 3] = _getvolume!(xex, f, SVector(x[m], x[m+1]),
                                            SVector(y[n], y[n+1]),
                                            SVector(coor[m, n, p-1][3], coor[m, n, p][3]))
    end

    w
end

# measure
function _getlength!(xex, f, x::SVector{2})
    t = SVector{2}(f(i) for i in x)

    val = x[2] - x[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return x[2] - ξ
    isnonnegative(t[2]) && return ξ - x[1]

    nothing
end

function _getarea!(xex, f, x::SVector{2}, y::SVector{2})
    t = SArray{Tuple{2,2}}(f(i, j) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    val * getcc(f, x0, h0, xex, Cint(2); nex=Cint.((0, 0)))
end

function _getvolume!(xex, f, x::SVector{2}, y::SVector{2}, z::SVector{2})
    t = SArray{Tuple{2,2,2}}(f(i, j, k) for i in x, j in y, k in z)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    val * getcc(f, x0, h0, xex, Cint(3); nex=Cint.((0, 0)))
end
