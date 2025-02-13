using CartesianGeometry
using Test
const T=Float64
using LinearAlgebra

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

@show Ws
@show length(Ws[1])

# Test ImplicitIntegration
levelset = (x) -> (x[1] - a)^2 + (x[2] - a)^2 - R^2

V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh, levelset)

@show W1
@show length(W1[1])

@test length(V) == length(V1)
@test all(@. isapprox(V, V1, atol=1e-7))
#@test cell_types1 == cell_types

#@test all(@. isapprox(Ws[1], W1[1], atol=1e-7))
#@test all(@. isapprox(Ws[2], W1[2], atol=1e-7))




cell_types = reshape(Ws[1], nx, ny)
cell_types1 = reshape(W1[1], nx, ny)

using CairoMakie

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = 1)
ax1 = Axis(fig[1, 2], aspect = 1)
hm = heatmap!(ax, cell_types, colormap = :thermal)
hm1 = heatmap!(ax1, cell_types1, colormap = :thermal)
Colorbar(fig[1, 3], hm)
Colorbar(fig[1, 4], hm1)

display(fig)
