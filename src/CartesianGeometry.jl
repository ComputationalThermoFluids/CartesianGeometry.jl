module CartesianGeometry

using StaticArrays
using Vofinit

const ArrayAbstract{N,T} = AbstractArray{T,N}

export HyperSphere
export mesh
export integrate

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
include("volume.jl")

end
