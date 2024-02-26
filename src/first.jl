droplast(this::Base.OneTo, n=1) = Base.OneTo(length(this)-n)
droplast(this::AbstractRange, n=1) = this[begin:end-n]

"""

    integrate(Tuple{0}, f, xyz, T=Float64)

Computes volume-specific (`Tuple{0}`) apertures of the first kind.

Returns a `Tuple` where

1. The first component is a `Vector{T}` that stores the wet volumes of each cell,
1. The second component is a `Vector{SVector{N,T}}` that stores the coordinates of the wet barycenters.

Wherever the moments can not be computed, the values are set to zero.

# Arguments

- `f`: the level set function.
- `xyz`: the Cartesian coordinates of the lattice nodes.

"""
function integrate(::Type{Tuple{0}}, f, xyz, T=Float64)
    N, len = length(xyz), prod(length.(xyz))

    v = Vector{T}(undef, len)
    bary = Vector{SVector{N,T}}(undef, len)

    integrate!((v, bary), Tuple{0}, f, xyz)
end

@generated function integrate!(moms, ::Type{Tuple{0}},
                               f, xyz::NTuple{N}) where {N}
    quote
        # axes
        input = only.(axes.(xyz))
        output = droplast.(input)

        # indices
        linear = LinearIndices(input)
        cartesian = CartesianIndices(output)

        (v, bary) = moms
        @ntuple($N, x) = xyz

        xex = zeros(Cdouble, 4)

        for index in cartesian
            n = linear[index]

            @ntuple($N, i) = Tuple(index)

            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))

            v[n] = @ncall($N, vofinit!, xex, f, y)
            bary[n] = @ncall($N, SVector, d -> xex[d])
        end

        # boundary conditions
        halo = EdgeIterator(CartesianIndices(input), cartesian)

        for index in halo
            n = linear[index]

            v[n] = zero(eltype(v))
            bary[n] = zero(eltype(bary))
        end

        return moms
    end
end

"""

    integrate(Tuple{1}, f, xyz, T=Float64)

Computes area-specific (`Tuple{1}`) apertures of the first kind.

Returns a `NTuple` where each element corresponds to  direction.

# Arguments

- `f`: the level set function.
- `xyz`: the Cartesian coordinates of the lattice nodes.

"""
function integrate(::Type{Tuple{1}}, f, xyz, T=Float64)
    len = prod(length.(xyz))

    moms = map(xyz) do _
        Vector{T}(undef, len)
    end

    integrate!(moms, Tuple{1}, f, xyz)
end

# 1D version
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{1})
    # axes
    input = only.(axes.(xyz))

    # indices
    linear = LinearIndices(input)

    x, = xyz

    xex = zeros(Cdouble, 4)

    # x faces
    output = input
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i,) = Tuple(index)
        moms[1][n] = vofinit!(xex, f, x[i])
    end

    return moms
end

# 2D version
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{2})
    # axes
    input = only.(axes.(xyz))

    # indices
    linear = LinearIndices(input)

    x, y = xyz

    xex = zeros(Cdouble, 4)

    # x faces
    output = input[1], droplast(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j) = Tuple(index)
        moms[1][n] = vofinit!(xex, f, x[i],
                                      SVector(y[j], y[j+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = zero(eltype(moms[1]))
    end

    # y faces
    output = droplast(input[1]), input[2]
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j) = Tuple(index)
        moms[2][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      y[j])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = zero(eltype(moms[2]))
    end

    return moms
end

# 3D version
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{3})
    # axes
    input = only.(axes.(xyz))

    # indices
    linear = LinearIndices(input)

    x, y, z = xyz

    xex = zeros(Cdouble, 4)

    # x faces
    output = input[1], droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j, k) = Tuple(index)
        moms[1][n] = vofinit!(xex, f, x[i],
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = zero(eltype(moms[1]))
    end

    # y faces
    output = droplast(input[1]), input[2], droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j, k) = Tuple(index)
        moms[2][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      y[j],
                                      SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = zero(eltype(moms[2]))
    end

    # z faces
    output = droplast(input[1]), droplast(input[2]), input[3]
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j, k) = Tuple(index)
        moms[3][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      z[k])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[3][n] = zero(eltype(moms[3]))
    end

    return moms
end
