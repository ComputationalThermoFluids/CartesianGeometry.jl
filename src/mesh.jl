function mesh(f, outer, inner)
    x = Vector{Float64}(undef, length(outer))

    lo = findfirst(isequal(first(inner)), outer)
    up = findfirst(isequal(last(inner)), outer)

    for i in eachindex(outer)
        x[i] = f((i-lo) / (up-lo))
    end

    x
end
