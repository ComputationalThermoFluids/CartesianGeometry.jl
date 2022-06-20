const ArrayAbstract{N,T} = AbstractArray{T,N}

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonnegative(x) = x ≥ zero(x)
isnonpositive(x) = x ≤ zero(x)

"""
    findin(subset, collection)

Return the indices of the elements of `subset` in `collection`. Return `nothing` if one or
more values is not in `collection`.

# Related PR and issues

- <https://github.com/JuliaLang/julia/pull/24673>
- <https://github.com/JuliaLang/julia/issues/30368>

"""
function findin(subset::AbstractUnitRange, collection::AbstractUnitRange)
    start = findfirst(isequal(first(subset)), collection)
    isnothing(start) && return nothing

    stop = findfirst(isequal(last(subset)), collection)
    isnothing(stop) && return nothing

    promote_type(typeof(subset), typeof(collection))(start, stop)
end
