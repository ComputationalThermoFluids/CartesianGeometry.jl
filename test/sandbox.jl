using CartesianGeometry

outer = (-1:10, -2:11)#, -1:18)
inner = (1:8, 1:8)#, 1:16)

xyz = mesh.(Ref(identity), outer, inner)

ranges = (; outer, inner)

levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

v = integrate(levelset, xyz, ranges)
