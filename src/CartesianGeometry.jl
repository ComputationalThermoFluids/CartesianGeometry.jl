module CartesianGeometry

using Base.Cartesian
using StaticArrays
using TiledIteration
using Vofinit

const ArrayAbstract{N,T} = AbstractArray{T,N}

export HyperSphere
export edge, center, mesh
export Dirichlet, apply!
export integrate

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
include("border.jl")
include("volume.jl")

end
