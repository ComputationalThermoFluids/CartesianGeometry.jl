using Statistics
using Plots
Plots.default(show = true)

# Function to compute wetted length along an edge
function compute_wetted_length(phi1, phi2, length)
    if phi1 >= 0 && phi2 >= 0
        return length
    elseif phi1 < 0 && phi2 < 0
        return 0.0
    else
        t = phi1 / (phi1 - phi2)
        if phi1 >= 0 && phi2 < 0
            return t * length
        else
            return (1 - t) * length
        end
    end
end

# Function to classify cells, find intersection points, and compute centroid and area
function marching_squares(x, y, phi::Array{Float64,2})
    nx = length(x) - 1  # Number of cells in x-direction
    ny = length(y) - 1  # Number of cells in y-direction

    cell_type = zeros(Int8, nx, ny)  # 0: empty, 1: full, -1: cut
    intersection_points = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}()
    centroids = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()
    V = zeros(Float64, nx, ny)
    interface_centroids = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()
    Ax = zeros(Float64, nx, ny)
    Ay = zeros(Float64, nx, ny)

    # Todo: Add Wx and Wy
    Wx = zeros(Float64, nx+1, ny)
    Wy = zeros(Float64, nx, ny+1)

    # Todo: Add Bx and By
    Bx = zeros(Float64, nx, ny)
    By = zeros(Float64, nx, ny)


    for j in 2:ny+1
        for i in 2:nx+1
            # Level set values at the corners
            phi_SW = phi[i-1,j-1]
            phi_SE = phi[i,j-1]
            phi_NE = phi[i,j]
            phi_NW = phi[i-1,j]

            # Determine the sign at each corner
            SW = sign(phi_SW)
            SE = sign(phi_SE)
            NE = sign(phi_NE)
            NW = sign(phi_NW)

            # Map signs to binary values: negative (-1) => 0, positive (1) => 1
            SW_bin = SW == -1 ? 0 : 1
            SE_bin = SE == -1 ? 0 : 1
            NE_bin = NE == -1 ? 0 : 1
            NW_bin = NW == -1 ? 0 : 1

            # Compute the isovalue
            iso = SW_bin + 2 * SE_bin + 4 * NE_bin + 8 * NW_bin

            # Coordinates of cell corners
            x1 = x[i-1]
            x2 = x[i]
            y1 = y[j-1]
            y2 = y[j]

            Δx = x2 - x1
            Δy = y2 - y1

            if iso == 0
                cell_type[i-1,j-1] = 0  # Empty cell

                # Store the centroid for this cell
                Cx = 0.5 * (x1 + x2)
                Cy = 0.5 * (y1 + y2)
                centroids[(i-1, j-1)] = (Cx, Cy)
                Ax[i-1,j-1] = 0.0
                Ay[i-1,j-1] = 0.0
                Bx[i-1,j-1] = 0.0
                By[i-1,j-1] = 0.0

                Wx[i-1,j-1] = 0.0
                Wy[i-1,j-1] = 0.0


            elseif iso == 15
                cell_type[i-1,j-1] = 1  # Full cell

                # Compute centroid and area for full cell
                Cx = 0.5 * (x1 + x2)
                Cy = 0.5 * (y1 + y2)
                area = Δx * Δy

                # Store the centroid and area for this cell
                centroids[(i-1, j-1)] = (Cx, Cy)
                V[i-1, j-1] = area
                Ax[i-1,j-1] = Δy  # Length of the right edge
                Ay[i-1,j-1] = Δx  # Length of the top edge
                Bx[i-1,j-1] = Δy  # Length of the horizontal center edge
                By[i-1,j-1] = Δx  # Length of the vertical center edge

                Wx[i-1,j-1] = Δy*Δx
                Wy[i-1,j-1] = Δy*Δx

            else
                cell_type[i-1,j-1] = -1  # Cut cell

                # Compute Ax and Ay for cut cells
                # Ax corresponds to the wetted length along the right edge
                Ax[i-1,j-1] = compute_wetted_length(phi_SE, phi_NE, Δy)

                # Ay corresponds to the wetted length along the top edge
                Ay[i-1,j-1] = compute_wetted_length(phi_NW, phi_NE, Δx)

                # Proceed with existing code for cut cells
                # Find intersection points along edges
                points = []

                # Edge connections
                edges = [
                    (phi_SW, phi_SE, (x1, y1), (x2, y1)),  # Edge 0: Bottom
                    (phi_SE, phi_NE, (x2, y1), (x2, y2)),  # Edge 1: Right
                    (phi_NE, phi_NW, (x2, y2), (x1, y2)),  # Edge 2: Top
                    (phi_NW, phi_SW, (x1, y2), (x1, y1))   # Edge 3: Left
                ]

                # Determine which corners are inside the interface
                inside_corners = []
                if phi_SW >= 0
                    push!(inside_corners, (x1, y1))
                end
                if phi_SE >= 0
                    push!(inside_corners, (x2, y1))
                end
                if phi_NE >= 0
                    push!(inside_corners, (x2, y2))
                end
                if phi_NW >= 0
                    push!(inside_corners, (x1, y2))
                end

                # Collect intersection points
                for edge in edges
                    phi1, phi2, (x1_edge, y1_edge), (x2_edge, y2_edge) = edge
                    if (phi1 >= 0 && phi2 < 0) || (phi1 < 0 && phi2 >= 0)
                        # There is an intersection along this edge
                        t = phi1 / (phi1 - phi2)
                        x_int = x1_edge + t * (x2_edge - x1_edge)
                        y_int = y1_edge + t * (y2_edge - y1_edge)
                        push!(points, (x_int, y_int))
                    end
                end

                # Compute the interface centroid (midpoint of the intersection segment)
                if length(points) == 2
                    (x_int1, y_int1) = points[1]
                    (x_int2, y_int2) = points[2]
                    Cx_interface = (x_int1 + x_int2) / 2
                    Cy_interface = (y_int1 + y_int2) / 2
                    interface_centroids[(i-1, j-1)] = (Cx_interface, Cy_interface)
                else
                    # Handle special cases where there are more than two intersection points
                    interface_centroids[(i-1, j-1)] = (NaN, NaN)  # Undefined for now
                end

                # Construct the polygon by combining inside corners and intersection points
                polygon_points = vcat(inside_corners, points)
                centroid_x = mean(first.(polygon_points))
                centroid_y = mean(last.(polygon_points))

                # Sort the points to form a proper polygon
                function angle(p)
                    return atan(p[2] - centroid_y, p[1] - centroid_x)
                end
                sorted_points = sort(polygon_points, by=angle)

                # Compute area and centroid using the shoelace formula
                n = length(sorted_points)
                area = 0.0
                Cx = 0.0
                Cy = 0.0

                for k in 1:n
                    xk, yk = sorted_points[k]
                    xk1, yk1 = sorted_points[mod1(k+1, n)]
                    cross = xk * yk1 - xk1 * yk
                    area += cross
                    Cx += (xk + xk1) * cross
                    Cy += (yk + yk1) * cross
                end
                area *= 0.5
                if area != 0.0
                    Cx /= (6 * area)
                    Cy /= (6 * area)
                else
                    Cx, Cy = centroid_x, centroid_y  # Default to mean if area is zero
                end

                # Store the intersection points, centroid, and area for this cell
                intersection_points[(i-1, j-1)] = sorted_points
                centroids[(i-1, j-1)] = (Cx, Cy)
                V[i-1, j-1] = abs(area)

                # Compute Wx and Wy for cut cells


                # Compute Bx and By for cut cells
                # Bx corresponds to the wetted length along the horizontal center edge
                # By corresponds to the wetted length along the vertical center edge
                # It's a line passing through the centroid of the cell 
                
            end
        end
    end
    return cell_type, intersection_points, centroids, V, interface_centroids, Ax, Ay, Wx, Wy, Bx, By
end

# Example usage
# Generate a mesh and a signed distance function for a circle
nx, ny = 10, 10  # Number of grid points
x = range(-1.0, 1.0, length=nx+1)
y = range(-1.0, 1.0, length=ny+1)
phi = zeros(Float64, nx+1, ny+1)

# Level set function for a circle centered at (0,0) with radius 0.5
for j in 1:ny+1
    for i in 1:nx+1
        phi[i,j] = (sqrt(x[i]^2 + y[j]^2) - 0.5)
    end
end

# Classify the cells, find intersection points, and compute centroids and V
@time cell_types, intersection_points, centroids, V, interface_centroids, Ax, Ay, Wx, Wy, Bx, By = marching_squares(x, y, phi)

# Display the results
println("Number of cells: ($(size(cell_types,1)), $(size(cell_types,2)))")
println("Number of centroids: $(length(centroids))")

# Visualize the cell types
heatmap(cell_types', aspect_ratio=:equal, c=:viridis, xlabel="x", ylabel="y", title="Cell Types")
readline()

# Visualize the cell V
heatmap(V', aspect_ratio=:equal, c=:viridis, xlabel="x", ylabel="y", title="Cell V")
readline()

# Visualize face capacities Ax and Ay
heatmap(Ax', aspect_ratio=:equal, c=:viridis, xlabel="x", ylabel="y", title="Face Capacity Ax")
readline()
heatmap(Ay', aspect_ratio=:equal, c=:viridis, xlabel="x", ylabel="y", title="Face Capacity Ay")
readline()

# Visualize the cell Wx and Wy
heatmap(Wx', aspect_ratio=:equal, c=:viridis, xlabel="x", ylabel="y", title="Cell Wx")
readline()
heatmap(Wy', aspect_ratio=:equal, c=:viridis, xlabel="x", ylabel="y", title="Cell Wy")
readline()


# Visualization of the classification and interface centroids
function plot_classification(x, y, cell_types, intersection_points, centroids, interface_centroids)
    nx, ny = size(cell_types)
    # Plot the grid nodes
    for xi in x
        for yj in y
            scatter!([xi], [yj], marker=:circle, color=:black)
        end
    end

    # Overlay the intersection points and centroids
    for ((i, j), points) in intersection_points
        # Plot the polygon
        x_coords = first.(points)
        y_coords = last.(points)
        plot!([x_coords; x_coords[1]], [y_coords; y_coords[1]], color=:red)
    end

    # Plot centroids
    for ((i, j), (Cx, Cy)) in centroids
        scatter!([Cx], [Cy], marker=:circle, color=:blue)
    end

    # Plot interface centroids
    for ((i, j), (Cx_int, Cy_int)) in interface_centroids
        if !isnan(Cx_int)
            scatter!([Cx_int], [Cy_int], marker=:diamond, color=:green)
        end
    end

    display(current())
end

#plot_classification(x, y, cell_types, intersection_points, centroids, interface_centroids)
#readline()
