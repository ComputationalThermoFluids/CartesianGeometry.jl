function integrate!(mom, ::Type{Tuple{P}}, f, xyz, domain) where {P}
    reshaped = reshape(mom)
    xyz = reshape.(xyz)

    _integrate!(reshaped, Tuple{P}, f, xyz, domain)

    mom
end

@generated function _integrate!(mom::ArrayAbstract{N}, ::Type{Tuple{0}},
                                f, xyz, domain) where {N}
    quote
        @ntuple($N, x) = xyz

        indices = CartesianIndices(domain)

        xex = zeros(Cdouble, 4)

        for index in indices
            @ntuple($N, i) = Tuple(index)

            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))

            vol = @ncall($N, vofinit!, xex, f, y)
            @nref($N, mom, i) = @ncall($N, SVector, vol, d -> xex[d])
        end

        mom
    end
end

function _integrate!(a::ArrayAbstract{2}, ::Type{Tuple{1}},
                     f, xyz, domains)
    (x,) = xyz

    nex = Cint.((0, 0))
    xex = zeros(Cdouble, 4)

    indices = CartesianIndices(domains[1])

    for index in indices
        (i,) = Tuple(index)
        a[i, 1] = vofinit!(xex, f, x[i])
    end

    a
end

function _integrate!(a::ArrayAbstract{3}, ::Type{Tuple{1}},
                     f, xyz, domains)
    (x, y) = xyz

    xex = zeros(Cdouble, 4)

    indices = CartesianIndices(domains[1])

    for index in indices
        (i, j) = Tuple(index)
        a[i, j, 1] = vofinit!(xex, f, x[i], SVector(y[j], y[j+1]))
    end

    indices = CartesianIndices(domains[2])

    for index in indices
        (i, j) = Tuple(index)
        a[i, j, 2] = vofinit!(xex, f, SVector(x[i], x[i+1]), y[j])
    end

    a
end

function _integrate!(a::ArrayAbstract{4}, ::Type{Tuple{1}},
                     f, xyz, domains)
    (x, y, z) = xyz

    xex = zeros(Cdouble, 4)

    indices = CartesianIndices(domains[1])

    for index in indices
        (i, j, k) = Tuple(index)
        a[i, j, k, 1] = vofinit!(xex, f, x[i], SVector(y[j], y[j+1]),
                                               SVector(z[k], z[k+1]))
    end

    indices = CartesianIndices(domains[2])

    for index in indices
        (i, j, k) = Tuple(index)
        a[i, j, k, 2] = vofinit!(xex, f, SVector(x[i], x[i+1]), y[j],
                                         SVector(z[k], z[k+1]))
    end

    indices = CartesianIndices(domains[3])

    for index in indices
        (i, j, k) = Tuple(index)
        a[i, j, k, 3] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                         SVector(y[j], y[j+1]), z[k])
    end

    a
end
