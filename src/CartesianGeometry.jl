module CartesianGeometry

using Base.Cartesian
using StaticArrays
using TiledIteration
using Vofinit

const ArrayAbstract{N,T} = AbstractArray{T,N}

export HyperSphere
export mesh
export Dirichlet, apply!
export integrate
export getvolume
export getsurface

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
include("border.jl")
include("volume.jl")
include("volume2.jl")
include("surface.jl")

end
