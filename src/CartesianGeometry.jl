module CartesianGeometry

using Base.Cartesian
using StaticArrays
using TiledIteration
using Vofinit
using CartesianArrays

#import Base: front, Fix1
#import Base: getindex

export HyperSphere
#export mesh
export collocated, staggered
#export Dirichlet, apply!
export integrate!#integrate
#export getvolume
#export getsurface

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
#include("border.jl")
include("volume.jl")
#include("volume2.jl")
include("surface.jl")
#include("surface2.jl")

end
