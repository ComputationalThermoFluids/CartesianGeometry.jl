using CartesianGeometry
using Test
const T=Float64

universe = (-1:11,)
node = (1:9,)

# define mesh
xyz = collocated.(identity, universe, node)

@assert all(@. isequal(length(xyz), length(universe)))

# define level set
const R = 0.25
const a = 0.5

levelset = HyperSphere(R, (a,))
#levelset = (x, _=0) -> (x - a)

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

@show cell_types
@show interface_length