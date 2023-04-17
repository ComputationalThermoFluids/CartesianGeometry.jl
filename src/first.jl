droplast(this::Base.OneTo, n=1) = Base.OneTo(length(this)-n)
droplast(this::AbstractRange, n=1) = this[begin:end-n]

"""

    integrate(Tuple{0}, f, absc, dom)

Computes volume-specific (`Tuple{0}`) apertures of the first kind.

Returns a `Tuple` which stores

1. The axes over which the apertures were computed.
1. The `Vector` that stores the apertures (length consistent with `absc`).

# Arguments

- `f`: the level set function.
- `absc`: the abscissas of the lattice nodes.

"""
function integrate(::Type{Tuple{0}}, f, absc, T=Float64)
    N, n = length(absc), prod(length.(absc))
    v = Vector{T}(undef, n)
    bary = Vector{SVector{N,T}}(undef, n)
    integrate!((v, bary), Tuple{0}, f, absc)
end

@generated function integrate!(moms, ::Type{Tuple{0}},
                               f, absc::NTuple{N}) where {N}
    quote
        # axes
        input = only.(axes.(absc))
        output = droplast.(input)

        # indices
        linear = LinearIndices(input)
        cartesian = CartesianIndices(output)

        (v, bary) = moms
        @ntuple($N, x) = absc

        xex = zeros(Cdouble, 4)

        for index in cartesian
            n = linear[index]

            @ntuple($N, i) = Tuple(index)

            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))

            v[n] = @ncall($N, vofinit!, xex, f, y)
            bary[n] = @ncall($N, SVector, d -> xex[d])
        end

        return output, moms...
    end
end

"""

    integrate(Tuple{1}, f, absc, dom)

Computes area-specific (`Tuple{1}`) apertures of the first kind.

"""
function integrate(::Type{Tuple{1}}, f, absc, T=Float64)
    avs = map(absc) do _
        Vector{T}(undef, prod(length.(absc)))
    end
    integrate!(avs, Tuple{1}, f, absc)
end

function integrate!(avs, ::Type{Tuple{1}}, f, absc::NTuple{1})
    # axes
    input = only.(axes.(absc))

    # indices
    linear = LinearIndices(input)

    (x,) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = input
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i,) = Tuple(index)
        avs[1][n] = vofinit!(xex, f, x[i])
    end

    axs = (output,)

    return axs, avs
end

function integrate!(avs, ::Type{Tuple{1}}, f, absc::NTuple{2})
    # axes
    input = only.(axes.(absc))

    # indices
    linear = LinearIndices(input)

    (x, y) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = input[1], droplast(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j) = Tuple(index)
        avs[1][n] = vofinit!(xex, f, x[i],
                                     SVector(y[j], y[j+1]))
    end

    axs = (output,)

    # y faces
    output = droplast(input[1]), input[2]
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j) = Tuple(index)
        avs[2][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                     y[j])
    end

    axs = (axs..., output)

    return axs, avs
end

function integrate!(avs, ::Type{Tuple{1}}, f, absc::NTuple{3})
    # axes
    input = only.(axes.(absc))

    # indices
    linear = LinearIndices(input)

    (x, y, z) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = input[1], droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j, k) = Tuple(index)
        avs[1][n] = vofinit!(xex, f, x[i],
                                     SVector(y[j], y[j+1]),
                                     SVector(z[k], z[k+1]))
    end

    axs = (output,)

    # y faces
    output = droplast(input[1]), input[2], droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j, k) = Tuple(index)
        avs[2][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                     y[j],
                                     SVector(z[k], z[k+1]))
    end

    axs = (axs..., output)

    # z faces
    output = droplast(input[1]), droplast(input[2]), input[3]
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        (i, j, k) = Tuple(index)
        avs[3][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                     SVector(y[j], y[j+1]),
                                     z[k])
    end

    axs = (axs..., output)

    return axs, avs
end
