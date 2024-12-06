using ImplicitIntegration
using Plots

# Define the mesh
nx, ny = 20, 20
dx, dy = 1/nx, 1/ny
x0, y0 = 0.0, 0.0

x_coords = [x0 + i*dx for i in 0:nx]
y_coords = [y0 + j*dy for j in 0:ny]

x_center = [0.5*(x_coords[i] + x_coords[i+1]) for i in 1:nx]
y_center = [0.5*(y_coords[j] + y_coords[j+1]) for j in 1:ny]

# Level-set function (example)
ϕ = (x) -> sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.25  # Example: unit circle

function compute_volume(ϕ, a, b)
    return ImplicitIntegration.integrate((x) -> 1, ϕ, a, b).val
end

# Compute the volume for each cell
volumes = zeros(nx, ny)
@time for i in 1:nx
    for j in 1:ny
        volumes[i, j] = compute_volume(ϕ, (x_coords[i], y_coords[j]), (x_coords[i+1], y_coords[j+1]))
    end
end

# Plot volumes
heatmap(x_center, y_center, volumes, aspect_ratio=1, c=:viridis, legend=false)

# Classify cells
cell_types = zeros(Int, nx, ny)

@time for i in 1:nx
    for j in 1:ny
        vol = volumes[i, j]
        if isapprox(vol, 0.0; atol=1e-8)
            cell_types[i, j] = 0  # Empty cell
        elseif isapprox(vol, dx * dy; atol=1e-8)
            cell_types[i, j] = 1  # Full cell
        else
            cell_types[i, j] = -1  # Cut cell
        end
    end
end

function compute_cell_centroid(ϕ, a, b)
    area = ImplicitIntegration.integrate((x) -> 1, ϕ, a, b).val
    x_centroid = ImplicitIntegration.integrate((x) -> x[1], ϕ, a, b).val / area
    y_centroid = ImplicitIntegration.integrate((x) -> x[2], ϕ, a, b).val / area
    return (x_centroid, y_centroid)
end

# Compute the centroid for each cell
centroids = Array{Tuple{Float64, Float64}, 2}(undef, nx, ny)

result = 0.0
@time for i in 1:nx
    for j in 1:ny
        a = (x_coords[i], y_coords[j])
        b = (x_coords[i+1], y_coords[j+1])
        result = compute_cell_centroid(ϕ, a, b)
        if isnan(result[1]) || isnan(result[2])
            centroids[i, j] = (0.5*(x_coords[i] + x_coords[i+1]), 0.5*(y_coords[j] + y_coords[j+1]))
        else
            centroids[i, j] = result
        end
    end
end

# Extract centroids for plotting
centroid_x = [c[1] for c in centroids]
centroid_y = [c[2] for c in centroids]

# Plot centroids
scatter!(repeat(x_center, 1, ny), repeat(y_center', nx, 1), markersize=5, c=:black, legend=false)
scatter!(reshape([c[1] for c in centroids], nx, ny), reshape([c[2] for c in centroids], nx, ny), markersize=5, c=:red, legend=false)

function compute_interface_centroid(ϕ, a, b)
    area = ImplicitIntegration.integrate((x) -> 1, ϕ, a, b,surface=true).val
    x_centroid = ImplicitIntegration.integrate((x) -> x[1], ϕ, a, b, surface=true).val / area
    y_centroid = ImplicitIntegration.integrate((x) -> x[2], ϕ, a, b, surface=true).val / area
    return (x_centroid, y_centroid)
end

# Compute the interface centroid for each cell
interface_centroids = Array{Union{Nothing, Tuple{Float64, Float64}}, 2}(undef, nx, ny)

@time for i in 1:nx
    for j in 1:ny
        a = (x_coords[i], y_coords[j])
        b = (x_coords[i+1], y_coords[j+1])
        interface_centroids[i, j] = compute_interface_centroid(ϕ, a, b)
    end
end

# Extract centroids for plotting
interface_centroid_x = [c[1] for c in interface_centroids if c !== nothing]
interface_centroid_y = [c[2] for c in interface_centroids if c !== nothing]

# Plot interface centroids
scatter!(interface_centroid_x, interface_centroid_y, markersize=5, c=:blue, legend=false)

# Function to compute staggered volume W_x at (i+1/2, j)
function compute_Wx(ϕ, xi, xip1, yj, yjp1)
    a = (xi, yj)
    b = (xip1, yjp1)
    return ImplicitIntegration.integrate((x) -> 1, ϕ, a, b).val
end

# Function to compute staggered volume W_y at (i, j+1/2)
function compute_Wy(ϕ, xi, xip1, yj, yjp1)
    a = (xi, yj)
    b = (xip1, yjp1)
    return ImplicitIntegration.integrate((x) -> 1, ϕ, a, b).val
end

# Compute staggered volumes W_x and W_y
Wx = zeros(nx+1, ny)
Wy = zeros(nx, ny+1)

# Compute W_x at x_{i+1/2}, y_j
for i in 1:nx+1
    xi = x_coords[max(i-1, 1)]
    xip1 = x_coords[min(i, nx+1)]
    for j in 1:ny
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        Wx[i, j] = compute_Wx(ϕ, xi, xip1, yj, yjp1)
    end
end

# Compute W_y at x_i, y_{j+1/2}
for i in 1:nx
    xi = x_coords[i]
    xip1 = x_coords[i+1]
    for j in 1:ny+1
        yj = y_coords[max(j-1, 1)]
        yjp1 = y_coords[min(j, ny+1)]
        Wy[i, j] = compute_Wy(ϕ, xi, xip1, yj, yjp1)
    end
end

# Function to compute Ax at x_{i+1/2}, y_j
function compute_Ax(ϕ, xi_half, yj_minus_half, yj_plus_half)
    # Define the fixed x coordinate
    x_fixed = xi_half
    # Define the lower and upper y coordinates
    y0 = yj_minus_half
    y1 = yj_plus_half
    # Perform 1D integration over y with x fixed
    a = (y0,)
    b = (y1,)
    ϕ_1d = (y) -> ϕ([x_fixed, y[1]])
    result = ImplicitIntegration.integrate((y) -> 1.0, ϕ_1d, a, b)
    return result.val
end

# Function to compute Ay at x_i, y_{j+1/2}
function compute_Ay(ϕ, xi_minus_half, xi_plus_half, yj_half)
    # Define the fixed y coordinate
    y_fixed = yj_half
    # Define the lower and upper x coordinates
    x0 = xi_minus_half
    x1 = xi_plus_half
    # Perform 1D integration over x with y fixed
    a = (x0,)
    b = (x1,)
    ϕ_1d = (x) -> ϕ([x[1], y_fixed])
    result = ImplicitIntegration.integrate((x) -> 1.0, ϕ_1d, a, b)
    return result.val
end

# Compute Ax and Ay for the mesh
Ax = zeros(nx+1, ny)
Ay = zeros(nx, ny+1)

# Compute Ax at x_{i+1/2}, y_j
for i in 1:nx+1
    xi_half = x_coords[i]
    for j in 1:ny
        yj_minus_half = y_coords[j]
        yj_plus_half = y_coords[j+1]
        Ax[i, j] = compute_Ax(ϕ, xi_half, yj_minus_half, yj_plus_half)
    end
end

# Compute Ay at x_i, y_{j+1/2}
for i in 1:nx
    xi_minus_half = x_coords[i]
    xi_plus_half = x_coords[i+1]
    for j in 1:ny+1
        yj_half = y_coords[j]
        Ay[i, j] = compute_Ay(ϕ, xi_minus_half, xi_plus_half, yj_half)
    end
end

heatmap(Ay, aspect_ratio=1, c=:viridis)

# Function to compute Bx at x_i, y_j
function compute_Bx(ϕ, xi_half, yj_minus_half, yj_plus_half)
    # Define the fixed x coordinate
    x_fixed = xi_half
    # Define the lower and upper y coordinates
    y0 = yj_minus_half
    y1 = yj_plus_half
    # Perform 1D integration over y with x fixed
    a = (y0,)
    b = (y1,)
    ϕ_1d = (y) -> ϕ([x_fixed, y[1]])
    result = ImplicitIntegration.integrate((y) -> 1.0, ϕ_1d, a, b)
    return result.val
end

# Function to compute By at x_i, y_j
function compute_By(ϕ, xi_minus_half, xi_plus_half, yj_half)
    # Define the fixed y coordinate
    y_fixed = yj_half
    # Define the lower and upper x coordinates
    x0 = xi_minus_half
    x1 = xi_plus_half
    # Perform 1D integration over x with y fixed
    a = (x0,)
    b = (x1,)
    ϕ_1d = (x) -> ϕ([x[1], y_fixed])
    result = ImplicitIntegration.integrate((x) -> 1.0, ϕ_1d, a, b)
    return result.val
end

# Compute Bx and By for the mesh
Bx = zeros(nx, ny)
By = zeros(nx, ny)

# Compute Bx at x_i, y_j
for i in 1:nx
    xi_half = centroid_x[i]
    for j in 1:ny
        yj_minus_half = y_coords[j]
        yj_plus_half = y_coords[j+1]
        Bx[i, j] = compute_Bx(ϕ, xi_half, yj_minus_half, yj_plus_half)
    end
end

# Compute By at x_i, y_j
for i in 1:nx
    xi_minus_half = x_coords[i]
    xi_plus_half = x_coords[i+1]
    for j in 1:ny
        yj_half = centroid_y[j]
        By[i, j] = compute_By(ϕ, xi_minus_half, xi_plus_half, yj_half)
    end
end

heatmap(Bx, aspect_ratio=1, c=:viridis)
