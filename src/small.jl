using StaticArrays

function find_small_cells(V::AbstractArray{T,N}, threshold::T) where {T,N}
    return CartesianIndices(V)[findall(v -> v > 0 && v < threshold, V)]
end

function find_neighbors(ind::CartesianIndex{N}, dims::NTuple{N,Int}) where N
    neighbor_indices = neighborsind(ind, dims)
    # Filtrer les voisins valides (dans les limites)
    return neighbor_indices
end

function neighborsind(ind::CartesianIndex{N}, dims::NTuple{N,Int}) where N
    neighbor_indices = CartesianIndex{N}[]
    for dim in 1:N
        for delta in (-1, 1)
            neighbor_coords = ntuple(i -> ind[i] + (i == dim ? delta : 0), N)
            if all(1 .≤ neighbor_coords .≤ dims)
                push!(neighbor_indices, CartesianIndex(neighbor_coords))
            end
        end
    end
    return neighbor_indices
end

function merge_cells!(V::AbstractArray{T,N}, bary::AbstractArray{SVector{N,T},N},
    idx::CartesianIndex{N}, neighbor_idx::CartesianIndex{N}) where {N,T}
V_new = V[idx] + V[neighbor_idx]
bary_new = (V[idx] * bary[idx] + V[neighbor_idx] * bary[neighbor_idx]) / V_new
V[neighbor_idx] = V_new
bary[neighbor_idx] = bary_new
V[idx] = 0.0
bary[idx] = zeros(SVector{N,T})
end
