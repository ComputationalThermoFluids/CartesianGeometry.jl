const ArrayAbstract{N,T} = AbstractArray{T,N}

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonnegative(x) = x ≥ zero(x)
isnonpositive(x) = x ≤ zero(x)
