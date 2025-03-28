using CartesianGeometry
using Test
const T=Float64
using LinearAlgebra

# define level set
const R = 0.21
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

@show As[1]
@show length(As[1])

# Test ImplicitIntegration
levelset = (x) -> (x[1] - a)^2 + (x[2] - a)^2 - R^2

V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh, levelset)

@show A1[1]
@show length(A1[1])

@test length(V) == length(V1)
@test all(@. isapprox(V, V1, atol=1e-7))
@test cell_types1 == cell_types
#@test all(@. isapprox(As[1], A1[1], atol=1e-7))

#@test all(@. isapprox(Ws[1], W1[1], atol=1e-7))
#@test all(@. isapprox(Ws[2], W1[2], atol=1e-7))

volumes = reshape(V, nx, ny)
volumes1 = reshape(V1, nx, ny)

using CairoMakie

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = 1)
ax1 = Axis(fig[1, 2], aspect = 1)
hm = heatmap!(ax, volumes, colormap = :thermal, label="V")
hm1 = heatmap!(ax1, volumes1, colormap = :thermal, label="V1")
Colorbar(fig[1, 3], hm)
Colorbar(fig[1, 4], hm1)
axislegend(ax, title = "V")
axislegend(ax1, title = "V1")

display(fig)


"""
cell_types = reshape(As[1], nx, ny)
cell_types1 = reshape(A1[1], nx, ny)

using CairoMakie

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = 1)
ax1 = Axis(fig[1, 2], aspect = 1)
hm = heatmap!(ax, cell_types, colormap = :thermal)
hm1 = heatmap!(ax1, cell_types1, colormap = :thermal)
Colorbar(fig[1, 3], hm)
Colorbar(fig[1, 4], hm1)

display(fig)

cell_centers_x = vec([c[1] for c in C_ω1])
cell_centers_y = vec([c[2] for c in C_ω1])

# some could be nothing if the cell is not an interface cell, handle it
interface_centers_x = vec([c[1] for c in C_γ1 if c !== nothing])
interface_centers_y = vec([c[2] for c in C_γ1 if c !== nothing])

# Prepare a figure to plot the cell centers.
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], title = "2D Cell Centers", xlabel = "x", ylabel = "y", aspect = DataAspect())
scatter!(ax, cell_centers_x, cell_centers_y, markersize = 10, color = :red)
scatter!(ax, interface_centers_x, interface_centers_y, markersize = 10, color = :blue)
fig
"""