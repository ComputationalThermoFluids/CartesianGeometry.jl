"""

    integrate(Tuple{0}, f, absc, dom)

Computes volume-specific (`Tuple{0}`) apertures of the first kind.

"""
function integrate(::Type{Tuple{0}}, f, absc, dom)
    T = SVector{length(absc)+1,Float64}
    mom = Vector{T}(undef, length(dom))
    integrate!(mom, Tuple{0}, f, absc, dom)
end

@generated function integrate!(mom, ::Type{Tuple{0}},
                               f, absc::NTuple{N}, dom) where {N}
    quote
        @ntuple($N, x) = absc

        xex = zeros(Cdouble, 4)

        for (i, index) in enumerate(dom)
            @ntuple($N, i) = Tuple(index)

            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))

            vol = @ncall($N, vofinit!, xex, f, y)
            mom[i] = @ncall($N, SVector, vol, d -> xex[d])
        end

        mom
    end
end

"""

    integrate(Tuple{1}, f, absc, dom)

Computes area-specific (`Tuple{1}`) apertures of the first kind.

"""
function integrate(::Type{Tuple{1}}, f, absc, doms)
    moms = map(doms) do dom
        Vector{Float64}(undef, length(dom))
    end
    integrate!(moms, Tuple{1}, f, absc, doms)
end

function integrate!(moms, ::Type{Tuple{1}}, f, absc::NTuple{1}, doms)
    (x,) = absc

    xex = zeros(Cdouble, 4)

    for (n, index) in enumerate(doms[1])
        (i,) = Tuple(index)
        moms[1][n] = vofinit!(xex, f, x[i])
    end

    moms
end

function integrate!(moms, ::Type{Tuple{1}}, f, absc::NTuple{2}, doms)
    (x, y) = absc

    xex = zeros(Cdouble, 4)

    for (n, index) in enumerate(doms[1])
        (i, j) = Tuple(index)
        moms[1][n] = vofinit!(xex, f, x[i],
                                      SVector(y[j], y[j+1]))
    end

    for (n, index) in enumerate(doms[2])
        (i, j) = Tuple(index)
        moms[2][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      y[j])
    end

    moms
end

function integrate!(moms, ::Type{Tuple{1}}, f, absc::NTuple{3}, doms)
    (x, y, z) = absc

    xex = zeros(Cdouble, 4)

    for (n, index) in enumerate(doms[1])
        (i, j, k) = Tuple(index)
        moms[1][n] = vofinit!(xex, f, x[i],
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]))
    end

    for (n, index) in enumerate(doms[2])
        (i, j, k) = Tuple(index)
        moms[2][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      y[j],
                                      SVector(z[k], z[k+1]))
    end

    for (n, index) in enumerate(doms[3])
        (i, j, k) = Tuple(index)
        moms[3][n] = vofinit!(xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      z[k])
    end

    moms
end
