function getsurface(f, xyz, ranges)
    (; outer) = ranges
    a = Vector{Float64}(undef, length(outer) * prod(length.(outer)))
    getsurface!(a, f, xyz, ranges)
end

function getsurface!(a, f, xyz, ranges)
    (; outer, stagger, center) = ranges

    ra = reshape(a, length.(outer)..., :)

    xyz = map(xyz, outer) do x, rng
        reshape(x, length(rng))
    end

    stagger = findin.(stagger, outer)
    center = findin.(center, outer)

    _getsurface!(ra, f, xyz, stagger, center)

    a
end

function _getsurface!(a::ArrayAbstract{2}, f, xyz, stagger, center)
    (x,) = xyz

    xex = zeros(Cdouble, 4)

    outer = CartesianIndices((stagger[1],))

    for mnp in outer
        (m,) = Tuple(mnp)

        a[m, 1] = _getpoint!(xex, f, x[m])
    end

    a
end

function _getsurface!(a::ArrayAbstract{3}, f, xyz, stagger, center)
    (x, y) = xyz

    xex = zeros(Cdouble, 4)

    # x-normal faces
    outer = CartesianIndices((stagger[1], center[2]))

    for mnp in outer
        (m, n) = Tuple(mnp)

        a[m, n, 1] = _getlength!(xex, f, x[m], SVector(y[n], y[n+1]))
    end

    # y-normal faces
    outer = CartesianIndices((center[1], stagger[2]))

    for mnp in outer
        (m, n) = Tuple(mnp)

        a[m, n, 2] = _getlength!(xex, f, SVector(x[m], x[m+1]), y[n])
    end

    a
end

function _getsurface!(a::ArrayAbstract{4}, f, xyz, stagger, center)
    (x, y, z) = xyz

    xex = zeros(Cdouble, 4)

    # x-normal faces
    outer = CartesianIndices((stagger[1], center[2], center[3]))

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        a[m, n, p, 1] = _getarea!(xex, f, x[m], SVector(y[n], y[n+1]),
                                                SVector(z[p], z[p+1]))
    end

    # y-normal faces
    outer = CartesianIndices((center[1], stagger[2], center[3]))

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        a[m, n, p, 2] = _getarea!(xex, f, SVector(x[m], x[m+1]), y[n],
                                          SVector(z[p], z[p+1]))
    end

    # z-normal faces
    outer = CartesianIndices((center[1], center[2], stagger[3]))

    for mnp in outer
        (m, n, p) = Tuple(mnp)

        a[m, n, p, 3] = _getarea!(xex, f, SVector(x[m], x[m+1]),
                                          SVector(y[n], y[n+1]), z[p])
    end

    a
end

function _getpoint!(xex, f, x)
    t = f(x)

    val = one(x)

    isnonpositive(t) && return val
    isnonnegative(t) && return zero(val)

    nothing
end

function _getlength!(xex, f, x::SVector{2}, y)
    t = SVector{2}(f(i, y) for i in x)

    val = x[2] - x[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return x[2] - ξ
    isnonnegative(t[2]) && return ξ - x[1]

    nothing
end

function _getlength!(xex, f, x, y::SVector{2})
    t = SVector{2}(f(x, j) for j in y)

    val = y[2] - y[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (y[2] * t[1] - y[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return y[2] - ξ
    isnonnegative(t[2]) && return ξ - y[1]

    nothing
end

function _getarea!(xex, f, x, y::SVector{2}, z::SVector{2})
    t = SMatrix{2,2}(f(x, j, k) for j in y, k in z)

    val = (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((y[1], z[1], zero(val)))
    h0 = Cdouble.((y[2]-y[1], z[2]-z[1], one(val)))

    val * getcc(x0, h0, xex, Cint(2); nex=Cint.((0, 0))) do y, z, _
        f(x, y, z)
    end
end

function _getarea!(xex, f, x::SVector{2}, y, z::SVector{2})
    t = SMatrix{2,2}(f(i, y, k) for i in x, k in z)

    val = (z[2] - z[1]) * (x[2] - x[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((z[1], x[1], zero(val)))
    h0 = Cdouble.((z[2]-z[1], x[2]-x[1], one(val)))

    val * getcc(x0, h0, xex, Cint(2); nex=Cint.((0, 0))) do z, x, _
        f(x, y, z)
    end
end

function _getarea!(xex, f, x::SVector{2}, y::SVector{2}, z)
    t = SMatrix{2,2}(f(i, j, z) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    val * getcc(x0, h0, xex, Cint(2); nex=Cint.((0, 0))) do x, y, _
        f(x, y, z)
    end
end
