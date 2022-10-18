# second kind capacities have second argument bary
function integrate!(w, ::Type{Tuple{P}}, f, xyz, bary, arg) where {P}
    reshaped = reshape(w)
    xyz = reshape.(xyz)
    bary = reshape(bary)

    _integrate!(reshaped, Tuple{P}, f, xyz, bary, arg)

    w
end

function _integrate!(w::ArrayAbstract{2}, ::Type{Tuple{0}}, f, xyz, bary, domains)
    (x,) = xyz

    nex = Cint.((0, 0))
    xex = zeros(Cdouble, 4)

    # x-normal faces
    indices = CartesianIndices(domains[1])

    for index in indices
        (i,) = Tuple(index)
        w[i, 1] = vofinit!(xex, f, SVector(bary[i-1][1], bary[i][1]); nex)
    end

    w
end

function _integrate!(w::ArrayAbstract{3}, ::Type{Tuple{0}}, f, xyz, bary, domains)
    (x, y) = xyz

    nex = Cint.((0, 0))
    xex = zeros(Cdouble, 4)

    # x-normal faces
    indices = CartesianIndices(domains[1])

    for index in indices
        (i, j) = Tuple(index)
        w[i, j, 1] = vofinit!(xex, f, SVector(bary[i-1, j][1], bary[i, j][1]),
                                      SVector(y[j], y[j+1]); nex)
    end

    # y-normal faces
    indices = CartesianIndices(domains[2])

    for index in indices
        (i, j) = Tuple(index)
        w[i, j, 2] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      SVector(bary[i, j-1][2], bary[i, j][2]); nex)
    end

    w
end

function _integrate!(w::ArrayAbstract{4}, ::Type{Tuple{0}}, f, xyz, bary, domains)
    (x, y, z) = xyz

    nex = Cint.((0, 0))
    xex = zeros(Cdouble, 4)

    # x-normal faces
    indices = CartesianIndices(domains[1])

    for index in indices
        (i, j, k) = Tuple(index)
        w[i, j, k, 1] = vofinit!(xex, f, SVector(bary[i-1, j, k][1], bary[i, j, k][1]),
                                         SVector(y[j], y[j+1]),
                                         SVector(z[k], z[k+1]); nex)
    end

    # y-normal faces
    indices = CartesianIndices(domains[2])

    for index in indices
        (i, j, k) = Tuple(index)
        w[i, j, k, 2] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                         SVector(bary[i, j-1, k][2], bary[i, j, k][2]),
                                         SVector(z[k], z[k+1]); nex)
    end

    # z-normal faces
    indices = CartesianIndices(domains[3])

    for index in indices
        (i, j, k) = Tuple(index)
        w[i, j, k, 3] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                         SVector(y[j], y[j+1]),
                                         SVector(bary[i, j, k-1][3], bary[i, j, k][3]); nex)
    end

    w
end

function _integrate!(b::ArrayAbstract{2}, ::Type{Tuple{1}}, f, xyz, bary, domain)
    (x,) = xyz

    xex = zeros(Cdouble, 4)

    indices = CartesianIndices(domain)

    # x-normal faces
    for index in indices
        (i,) = Tuple(index)
        b[i, 1] = vofinit!(xex, f, bary[i][1])
    end

    b
end

function _integrate!(b::ArrayAbstract{3}, ::Type{Tuple{1}}, f, xyz, bary, domain)
    (x, y) = xyz

    xex = zeros(Cdouble, 4)

    indices = CartesianIndices(domain)

    # x-normal faces
    for index in indices
        (i, j) = Tuple(index)
        b[i, j, 1] = vofinit!(xex, f, bary[i, j][1], SVector(y[j], y[j+1]))
    end

    # y-normal faces
    for index in indices
        (i, j) = Tuple(index)
        b[i, j, 2] = vofinit!(xex, f, SVector(x[i], x[i+1]), bary[i, j][2])
    end

    b
end

function _integrate!(b::ArrayAbstract{4}, ::Type{Tuple{1}}, f, xyz, bary, domain)
    (x, y, z) = xyz

    xex = zeros(Cdouble, 4)

    indices = CartesianIndices(domain)

    # x-normal faces
    for index in indices
        (i, j, k) = Tuple(index)
        b[i, j, k, 1] = vofinit!(xex, f, bary[i, j, k][1], SVector(y[j], y[j+1]),
                                                           SVector(z[k], z[k+1]))
    end

    # y-normal faces
    for index in indices
        (i, j, k) = Tuple(index)
        b[i, j, k, 2] = vofinit!(xex, f, SVector(x[i], x[i+1]), bary[i, j, k][2],
                                         SVector(z[k], z[k+1]))
    end

    # z-normal faces
    for index in indices
        (i, j, k) = Tuple(index)
        b[i, j, k, 3] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                         SVector(y[j], y[j+1]), bary[i, j, k][3])
    end

    b
end
