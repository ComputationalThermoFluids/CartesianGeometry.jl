using CartesianGeometry

outer = (-1:10,)
inner = (1:8,)

xyz = mesh.(Ref(identity), outer, inner)

ranges = (; outer, inner)

v = integrate(xyz, ranges) do x
    x - one(x) / 2
end

