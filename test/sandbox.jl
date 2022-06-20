using CartesianArrays
using CartesianGeometry
using StaticArrays

import Base: OneTo

const T = Float64

# "infinite" range for mesh abscissas
# to not think about it
infinite = (-1:11, -1:19)

# primary ranges
outer = (0:10, 0:18)
node = (1:9, 1:17)
center = (1:8, 1:16)

# secondary ranges
faces = ntuple(length(outer)) do i
    ntuple(j -> j==i ? node[i] : center[j], length(outer))
end
outers = ntuple(length(outer)) do _
    outer
end

xyz = map(infinite, node) do o, i
    a = staggered(identity, o, i)
    CartesianArray(a, tuple(o))
end

xyz_ = map(infinite, center) do o, i
    a = centered(identity, o, i)
    CartesianArray(a, tuple(o))
end

# geometry
levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

v = CartesianArray{SVector{length(outer)+1,T}}(undef, outer)
integrate!(v, Tuple{0}, levelset, xyz, outer)

a = CartesianArray{T}(undef, (outer..., OneTo(length(outer))))
integrate!(a, Tuple{1}, levelset, xyz, outers)

#b = CartesianArray{T}(undef, (outer..., OneTo(length(outer))))
#integrate!(b, Tuple{1}, levelset, xyz, v, outer)

#= from here
(v, coor) = integrate(Tuple{0,1}, levelset, sxyz, infinite, couter, cinner)
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
