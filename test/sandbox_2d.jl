using CartesianGeometry
using Test
const T=Float64

# define level set
const R = 0.2
const a = 0.5

nx, ny = 10, 10
dx, dy = 1.0/nx, 1.0/ny
x0, y0 = 0.0, 0.0

mesh = ([x0 + i*dx for i in 0:nx-1], [y0 + i*dy for i in 0:ny-1])

levelset = (x, y, _=0) -> (x - a)^2 + (y - a)^2 - R^2

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, mesh, T, zero)
As = integrate(Tuple{1}, levelset, mesh, T, zero)

Ws = integrate(Tuple{0}, levelset, mesh, T, zero, bary)
Bs = integrate(Tuple{1}, levelset, mesh, T, zero, bary)

@show V

# Test ImplicitIntegration
levelset = (x) -> (x[1] - a)^2 + (x[2] - a)^2 - R^2

V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh, levelset)

@show V1

#@test cell_types1 == cell_types
