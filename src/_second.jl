dropfirst(this::AbstractRange, n=1) = this[begin+n:end]

"""

    integrate(T, f, absc, axs, bary)

Computes volume- (`T=Tuple{0}`) and surface-specific (`T=Tuple{1}`) apertures of the second kind.

# Arguments

- `axs`: the indices of `absc` for which `bary` is defined.

"""
function integrate(T::Type{<:Tuple},
                   f, absc, axs, bary, S=eltype(eltype(bary)))
    avs = map(absc) do _
        similar(bary, S)
    end
    integrate!(avs, T, f, absc, axs, bary)
end

function integrate!(avs, ::Type{Tuple{0}},
                    f, absc::NTuple{1}, axs, bary)
    input = only.(axes.(absc))
    linear = LinearIndices(input)
    cartesian = CartesianIndices(input)

    (x,) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = (dropfirst(axs[1]),)
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i,) = Tuple(index)
        left = i-1
        right = n
        avs[1][n] = vofinit!(xex, f,
                             SVector(bary[left][1], bary[right][1]))
    end

    res = (output,)

    return res, avs
end

function integrate!(avs, ::Type{Tuple{0}},
                    f, absc::NTuple{2}, axs, bary)
    input = only.(axes.(absc))
    linear = LinearIndices(input)
    cartesian = CartesianIndices(input)

    (x, y) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = dropfirst(axs[1]), droplast(input[2])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j) = Tuple(index)
        left = linear[CartesianIndex(i-1, j)]
        right = n
        avs[1][n] = vofinit!(xex, f,
                             SVector(bary[left][1], bary[right][1]),
                             SVector(y[j], y[j+1]))
    end

    res = (output,)

    # y faces
    output = droplast(input[1]), dropfirst(axs[2])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j) = Tuple(index)
        left = linear[CartesianIndex(i, j-1)]
        right = n
        avs[2][n] = vofinit!(xex, f,
                             SVector(x[i], x[i+1]),
                             SVector(bary[left][2], bary[right][2]))
    end

    res = (res..., output)

    return res, avs
end

function integrate!(avs, ::Type{Tuple{0}},
                    f, absc::NTuple{3}, axs, bary)
    input = only.(axes.(absc))
    linear = LinearIndices(input)
    cartesian = CartesianIndices(input)

    (x, y, z) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = dropfirst(axs[1]), droplast(input[2]), droplast(input[3])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j, k) = Tuple(index)
        left = linear[CartesianIndex(i-1, j, k)]
        right = n
        avs[1][n] = vofinit!(xex, f,
                              SVector(bary[left][1], bary[right][1]),
                              SVector(y[j], y[j+1]),
                              SVector(z[k], z[k+1]))
    end

    res = (output,)

    # y faces
    output = droplast(input[1]), dropfirst(axs[2]), droplast(input[3])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j, k) = Tuple(index)
        left = linear[CartesianIndex(i, j-1, k)]
        right = n
        avs[2][n] = vofinit!(xex, f,
                             SVector(x[i], x[i+1]),
                             SVector(bary[left][2], bary[right][2]),
                             SVector(z[k], z[k+1]))
    end

    res = (res..., output)

    # z faces
    output = droplast(input[1]), droplast(input[2]), dropfirst(axs[3])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j, k) = Tuple(index)
        left = linear[CartesianIndex(i, j, k-1)]
        right = n
        avs[3][n] = vofinit!(xex, f,
                             SVector(x[i], x[i+1]),
                             SVector(y[j], y[j+1]),
                             SVector(bary[left][3], bary[right][3]))
    end

    res = (res..., output)

    return res, avs
end

function integrate!(avs, ::Type{Tuple{1}},
                    f, absc::NTuple{1}, axs, bary)
    input = only.(axes.(absc))
    linear = LinearIndices(input)
    cartesian = CartesianIndices(input)

    (x,) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = (axs[1],)
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i,) = Tuple(index)
        avs[1][n] = vofinit!(xex, f, bary[i][1])
    end

    res = (output,)

    return res, avs
end

function integrate!(avs, ::Type{Tuple{1}},
                    f, absc::NTuple{2}, axs, bary)
    input = only.(axes.(absc))
    linear = LinearIndices(input)
    cartesian = CartesianIndices(input)

    (x, y) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = axs[1], droplast(input[2])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j) = Tuple(index)
        avs[1][n] = vofinit!(xex, f,
                             bary[n][1],
                             SVector(y[j], y[j+1]))
    end

    res = (output,)

    # y faces
    output = droplast(input[1]), axs[2]
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j) = Tuple(index)
        avs[2][n] = vofinit!(xex, f,
                             SVector(x[i], x[i+1]),
                             bary[n][2])
    end

    res = (res..., output)

    return res, avs
end

function integrate!(avs, ::Type{Tuple{1}},
                    f, absc::NTuple{3}, axs, bary)
    input = only.(axes.(absc))
    linear = LinearIndices(input)
    cartesian = CartesianIndices(input)

    (x, y, z) = absc

    xex = zeros(Cdouble, 4)

    # x faces
    output = axs[1], droplast(input[2]), droplast(input[3])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j, k) = Tuple(index)
        avs[1][n] = vofinit!(xex, f,
                             bary[n][1],
                             SVector(y[j], y[j+1]),
                             SVector(z[k], z[k+1]))
    end

    res = (output,)

    # y faces
    output = droplast(input[1]), axs[2], droplast(input[3])
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j, k) = Tuple(index)
        avs[2][n] = vofinit!(xex, f,
                             SVector(x[i], x[i+1]),
                             bary[n][2],
                             SVector(z[k], z[k+1]))
    end

    res = (res..., output)

    # z faces
    output = droplast(input[1]), droplast(input[2]), axs[3]
    indices = CartesianIndices(output)

    for index in indices
        n = linear[index]
        (i, j, k) = Tuple(index)
        avs[3][n] = vofinit!(xex, f,
                             SVector(x[i], x[i+1]),
                             SVector(y[j], y[j+1]),
                             bary[n][3])
    end

    res = (res..., output)

    return res, avs
end
