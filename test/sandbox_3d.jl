using CartesianGeometry
using Test
const T=Float64

# Test 3D
universe = (-1:11, -1:11, -1:11)
node = (1:9, 1:9, 1:9)

xyz = collocated.(identity, universe, node)

@assert all(@. isequal(length(xyz), length(universe)))

# define level set
const R = 0.25
const a, b, c = 0.5, 0.5, 0.5

levelset = HyperSphere(R, (a, b, c))

V, bary, interface_area = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

@test typeof(Bs) == Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

@test size(V) == size(bary) == size(interface_area)

@show interface_area