using CartesianGeometry

const T=Float64

xyz = (0:0.1:1, 0:0.1:1)

levelset = (x, y, _=0) -> x^2 + y^2 - 1

V, bary = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

