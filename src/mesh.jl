#const edge = Val{true}
#const center = Val{false}
#
# First and last points on boundaries
function mesh(::Type{Val{true}}, f, outer, inner)
    x = Vector{Float64}(undef, length(outer))

    lo = findfirst(isequal(first(inner)), outer)
    up = findfirst(isequal(last(inner)), outer)

    for i in eachindex(outer)
        x[i] = f((i-lo) / (up-lo))
    end

    x
end

# First and last points half a point away from boundaries
function mesh(::Type{Val{false}}, f, outer, inner)
    x = Vector{Float64}(undef, length(outer))

    lo = findfirst(isequal(first(inner)), outer)
    up = findfirst(isequal(last(inner)), outer)

    for i in eachindex(outer)
        x[i] = f((2(i-lo) + 1) / 2(up-lo+1))
    end

    x
end
