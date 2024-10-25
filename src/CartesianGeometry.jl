module CartesianGeometry

using Base.Cartesian
using StaticArrays
using TiledIteration
using Vofinit
#using CartesianCore
#using CartesianArrays

export nan
export HyperSphere
export collocated, staggered
export integrate, integrate!
export find_small_cells, find_neighbors, merge_cells!, neighborsind

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
include("vofinit.jl")
include("first.jl")
include("second.jl")
include("small.jl")

end
