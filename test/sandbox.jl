using StaticArrays
using CartesianArrays
using CartesianGeometry

import Base: OneTo

const T = Float64

#=== Domains ===#
universe = (-1:11, -1:19)

halos = ((-1:10, -1:18),
         (0:10, 0:18))
node = (1:9, 1:17)
center = (1:8, 1:16)

#=== Mesh ===#
xyz = collocated.(Ref(identity), universe, node)

#=== Geometry ===#
levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

#=== Capacities (1st kind) ===#
#--- Volume ---#
mom = CVector{SVector{length(halos[1])+1,T}}(undef, halos[1])
integrate!(mom, Tuple{0}, levelset, xyz, halos[1])

v = first.(mom)
bary = deleteat.(mom, 1)
mom = nothing

#--- Surface ---#
domains = ntuple(_ -> halos[1], length(halos[1]))
a = CVector{T}(undef, (halos[1]..., OneTo(length(halos[1]))))
integrate!(a, Tuple{1}, levelset, xyz, domains)

#=== Capacities (2nd kind) ===#
#--- Volume ---#
domains = ntuple(_ -> halos[2], length(halos[2]))
w = CVector{T}(undef, (halos[2]..., OneTo(length(halos[2]))))
integrate!(w, Tuple{0}, levelset, xyz, bary, domains)

#--- Surface ---#
b = CVector{T}(undef, (halos[1]..., OneTo(length(halos[1]))))
integrate!(b, Tuple{1}, levelset, xyz, bary, halos[1])
