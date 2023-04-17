module CartesianGeometry

using Base.Cartesian
using StaticArrays
using Vofinit
#using CartesianCore
#using CartesianArrays

#export HyperSphere
export collocated, staggered
export integrate, integrate!

include("utils.jl")
#include("zoo.jl")
include("mesh.jl")
include("vofinit.jl")
include("first.jl")
#include("second.jl")
include("_second.jl")

end
