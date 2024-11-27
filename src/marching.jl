# Function to classify cells, find intersection points, and compute centroid and area
function marching_squares(x, y, phi::Array{Float64,2})
    nx = length(x) - 1  # Number of cells in x-direction
    ny = length(y) - 1  # Number of cells in y-direction

    cell_type = zeros(Int8, nx, ny)  # 0: empty, 1: full, -1: cut
    intersection_points = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}()
    centroids = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()
    areas = Dict{Tuple{Int, Int}, Float64}()
    interface_centroids = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()

    for j in 2:ny
        for i in 2:nx
            # Bilinear interpolation to compute level set values at the corners
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

            # Classify the cell based on the isovalue
            if iso == 0
                cell_type[i-1,j-1] = 0  # Empty cell

                # Store the centroid for this cell
                Cx = 0.5 * (x[i-1] + x[i])
                Cy = 0.5 * (y[j-1] + y[j])
                centroids[(i-1, j-1)] = (Cx, Cy)
            elseif iso == 15
                cell_type[i-1,j-1] = 1  # Full cell

                # Compute centroid and area for full cell : Regular quadrilateral
                x1 = x[i-1]
                x2 = x[i]
                y1 = y[j-1]
                y2 = y[j]
                Cx = 0.5 * (x1 + x2)
                Cy = 0.5 * (y1 + y2)
                area = (x2 - x1) * (y2 - y1)

                # Store the centroid and area for this cell
                centroids[(i-1, j-1)] = (Cx, Cy)

                areas[(i-1, j-1)] = area

            else
                cell_type[i-1,j-1] = -1  # Cut cell

                # Find intersection points along edges
                points = []

                # Coordinates of cell corners
                x1 = x[i-1]
                x2 = x[i]
                y1 = y[j-1]
                y2 = y[j]

                # Level set values at corners
                phi_corners = [phi_SW, phi_SE, phi_NE, phi_NW]

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
                        t = phi1 / (phi1 - phi2)  # Fraction along the edge
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
                # The order of points is important; we can sort them by angle from the centroid
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
                areas[(i-1, j-1)] = abs(area)
            end
        end
    end
    return cell_type, intersection_points, centroids, areas, interface_centroids
end
