function vofinit!(xex, f, x::Number)
    t = f(x)

    val = one(x)

    isnonpositive(t) && return val
    isnonnegative(t) && return zero(val)

    nothing
end

"""

1d implementation complies with Vofi 2.0's API (see 2d and 3d).

!!! note "Convention"

    Even if cell is empty, return full centroid coordinates.

"""
function vofinit!(xex, f, x::SVector; nex=Cint.((1, 1)))
    t = SVector{2}(f(i) for i in x)

    val = x[2] - x[1]

    if all(isnonpositive, t)
        isone(first(nex)) && (xex[1] = sum(x) / 2)
        return val
    end

    if all(isnonnegative, t)
        isone(first(nex)) && (xex[1] = sum(x) / 2)
        return zero(val)
    end

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    if isnonnegative(t[1])
        isone(first(nex)) && (xex[1] = (ξ + x[2]) / 2)
        return x[2] - ξ
    end

    if isnonnegative(t[2])
        isone(first(nex)) && (xex[1] = (x[1] + ξ) / 2)
        return ξ - x[1]
    end

    nothing
end

function vofinit!(xex, f, x::SVector{2}, y::Number)
    t = SVector{2}(f(i, y) for i in x)

    val = x[2] - x[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return x[2] - ξ
    isnonnegative(t[2]) && return ξ - x[1]

    nothing
end

function vofinit!(xex, f, x::Number, y::SVector{2})
    t = SVector{2}(f(x, j) for j in y)

    val = y[2] - y[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (y[2] * t[1] - y[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return y[2] - ξ
    isnonnegative(t[2]) && return ξ - y[1]

    nothing
end

"""

Call Vofi 2.0 for exact integration.

"""
function vofinit!(xex, f, x::SVector, y::SVector; nex=Cint.((1, 1)))
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

    val * getcc(f, x0, h0, xex, Cint(2); nex)
end

function vofinit!(xex, f, x::Number, y::SVector{2}, z::SVector{2})
    t = SMatrix{2,2}(f(x, j, k) for j in y, k in z)

    val = (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((y[1], z[1], zero(val)))
    h0 = Cdouble.((y[2]-y[1], z[2]-z[1], one(val)))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(2); nex) do y, z, _
        f(x, y, z)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::Number, z::SVector{2})
    t = SMatrix{2,2}(f(i, y, k) for i in x, k in z)

    val = (z[2] - z[1]) * (x[2] - x[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((z[1], x[1], zero(val)))
    h0 = Cdouble.((z[2]-z[1], x[2]-x[1], one(val)))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(2); nex) do z, x, _
        f(x, y, z)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::Number)
    t = SMatrix{2,2}(f(i, j, z) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(2); nex) do x, y, _
        f(x, y, z)
    end
end

"""

Call Vofi 2.0 for exact integration.

"""
function vofinit!(xex, f, x::SVector, y::SVector, z::SVector; nex=Cint.((1, 1)))
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

    val * getcc(f, x0, h0, xex, Cint(3); nex)
end
