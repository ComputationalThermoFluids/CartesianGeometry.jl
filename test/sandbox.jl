using StaticArrays
using Accessors
using CartesianArrays
using CartesianGeometry

import Base: OneTo

const T = Float64

# "universe" domain for mesh abscissas
# to not think about it
universe = (-1:11, -1:19)

# primary domains
outer = (0:10, 0:18)
node = (1:9, 1:17)
center = (1:8, 1:16)

# secondary domains
#faces = ntuple(length(outer)) do i
#    ntuple(j -> j==i ? node[i] : center[j], length(outer))
#end
faces = ntuple(length(outer)) do d
    @set center[d] = node[d]
end
outers = ntuple(length(outer)) do _
    outer
end

xyz = collocated.(Ref(identity), universe, node)
#xyz_ = staggered.(Ref(identity), universe, center)

# geometry
levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

v = CVector{SVector{length(outer)+1,T}}(undef, outer)
integrate!(v, Tuple{0}, levelset, xyz, outer)

a = CVector{T}(undef, (outer..., OneTo(length(outer))))
integrate!(a, Tuple{1}, levelset, xyz, outers)

#b = CartesianArray{T}(undef, (outer..., OneTo(length(outer))))
#integrate!(b, Tuple{1}, levelset, xyz, v, outer)

#= from here
(v, coor) = integrate(Tuple{0,1}, levelset, sxyz, universe, couter, cinner)
to there =#

#    mesh.(Ref(Val{true}), Ref(identity), ranges, stagger)
#=
outer = (-1:10, -2:11)#, -1:10)#, -1:18)
stagger = (1:9, 1:9)#, 1:9)#, 1:16)

xyz = mesh.(Ref(Val{true}), Ref(identity), outer, stagger)

center = (1:8, 1:8)#, 1:8)
ranges = (; outer, inner=center)

levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

v, coor = integrate(levelset, xyz, ranges)

xyz2 = mesh.(Ref(Val{false}), Ref(identity), outer, center)

homogeneous(::T...) where {T} = zero(T)

apply!(v, Dirichlet(homogeneous, xyz2, ranges))

apply!(coor, Dirichlet(SVector{length(xyz),eltype(eltype(coor))}, xyz2, ranges))

ranges = (; outer, stagger, center)
a = getsurface(levelset, xyz, ranges)
w = getvolume(levelset, xyz, coor, ranges)
b = getsurface(levelset, xyz, coor, (; outer, inner=center))

=#
