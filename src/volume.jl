#=
"""

- `xyz`: abscissas of staggered edges;
- `stagger`, `center`: outer ranges;
- `inner`: centered inner range.

"""
function integrate(::Type{Tuple{0,1}}, f, xyz, stagger, center, inner)
    d = length(center)
    n = prod(length.(center))

    mom = (Vector{Float64}(undef, n),
           Vector{SVector{d,Float64}}(undef, n))

    integrate!(mom, f, xyz, stagger, center, inner)
end

function integrate!(mom, f, xyz, stagger, center, inner)
    dims = length.(center)

    rmom = map(mom) do el
        reshape(el, dims...)
    end

    rxyz = map(xyz, stagger) do el, range
        reshape(el, length(range))
    end

    stagger = findin.(inner, stagger)
    center = findin.(inner, center)

    _integrate!(rmom, f, rxyz, stagger, center)

    mom
end

=#

function integrate!(mom, ::Type{Tuple{0}}, f, xyz, ranges=getranges(mom, 1))
    center = findin.(ranges, getranges(mom, 1))
    reshaped = reshape(mom)

    node = findin.(ranges, only.(getranges.(xyz, 1)))
    xyz = reshape.(xyz)

    _integrate!(reshaped, Tuple{0}, f, xyz, node, center)

    mom
end

@generated function _integrate!(mom::ArrayAbstract{N}, ::Type{Tuple{0}},
                                f, xyz, node, center) where {N}
    quote
        @ntuple($N, x) = xyz

        node = CartesianIndices(node)
        center = CartesianIndices(center)

        xex = zeros(Cdouble, 4)

        for (nind, cind) in zip(node, center)
            @ntuple($N, n) = Tuple(nind)
            @ntuple($N, c) = Tuple(cind)

            @nextract($N, y, d -> SVector(x_d[n_d], x_d[n_d+1]))

            vol = @ncall($N, _vofi_integrate!, xex, f, y)
            @nref($N, mom, c) = @ncall($N, SVector, vol, d -> xex[d])
        end

        mom
    end
end

#=
@generated function _integrate!(mom::NTuple{2,ArrayAbstract{N}},
                                f, xyz, stagger, center) where {N}
    quote
        (volume, centroid) = mom
        @ntuple($N, x) = xyz

        stagger = CartesianIndices(stagger)
        center = CartesianIndices(center)

        xex = zeros(Cdouble, 4)

        for (sind, cind) in zip(stagger, center)
            @ntuple($N, s) = Tuple(sind)
            @ntuple($N, c) = Tuple(cind)

            @nextract($N, y, d -> SVector(x_d[s_d], x_d[s_d+1]))

            @nref($N, volume, c) = @ncall($N, _vofi_integrate!, xex, f, y)
            @nref($N, centroid, c) = @ncall($N, SVector, d -> xex[d])
        end

        mom
    end
end

=#

"""

1d implementation complies with Vofi 2.0's API (see 2d and 3d).

!!! note "Convention"

    Even if cell is empty, return full centroid coordinates.

"""
function _vofi_integrate!(xex, f, x)
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

"""

Call Vofi 2.0 for exact integration.

"""
function _vofi_integrate!(xex, f, x, y)
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

"""

Call Vofi 2.0 for exact integration.

"""
function _vofi_integrate!(xex, f, x, y, z)
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
