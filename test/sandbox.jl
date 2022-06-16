using CartesianGeometry
using StaticArrays

outer = (-1:10, -2:11, -1:10)#, -1:18)
stagger = (1:9, 1:9, 1:9)#, 1:16)

xyz = mesh.(Ref(Val{true}), Ref(identity), outer, stagger)

center = (1:8, 1:8, 1:8)
ranges = (; outer, inner=center)

levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

v, coor = integrate(levelset, xyz, ranges)

xyz2 = mesh.(Ref(Val{false}), Ref(identity), outer, center)

homogeneous(::T...) where {T} = zero(T)

apply!(v, Dirichlet(homogeneous, xyz2, ranges))

apply!(coor, Dirichlet(SVector{length(xyz),eltype(eltype(coor))}, xyz2, ranges))

ranges = (; outer, stagger, center)
a = getsurface(levelset, xyz, ranges)
