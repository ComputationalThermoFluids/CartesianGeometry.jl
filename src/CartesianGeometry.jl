module CartesianGeometry

const ArrayAbstract{N,T} = AbstractArray{T,N}

export mesh
export integrate

include("utils.jl")
include("mesh.jl")
include("volume.jl")

end
