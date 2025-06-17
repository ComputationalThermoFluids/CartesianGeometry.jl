module SimpleVOF

export getcelltype, getcc
export compute_interface_length
export compute_interface_area

# 2D functions
"""
    getcelltype(ls, xcell, ycell)

Retourne le type de cellule pour une cellule délimitée par les coordonnées données:
- 1 : cellule pleine (level set < 0 partout)
- 0 : cellule vide (level set > 0 partout)
- -1 : cellule coupée (level set change de signe)
"""
function getcelltype(ls, xcell, ycell)
    # Évaluer la level set aux 4 sommets de la cellule
    φ = [
        ls(xcell[1], ycell[1]),  # coin bas-gauche
        ls(xcell[2], ycell[1]),  # coin bas-droite
        ls(xcell[2], ycell[2]),  # coin haut-droite
        ls(xcell[1], ycell[2])   # coin haut-gauche
    ]
    
    if all(φ .< 0)
        return 1  # cellule pleine
    elseif all(φ .> 0)
        return 0  # cellule vide
    else
        return -1  # cellule coupée
    end
end

"""
    find_zero(ls, p1, p2)

Finds the zero of the level set function on segment [p1, p2] using linear interpolation.
"""
function find_zero(ls, p1, p2)
    x1, y1 = p1
    x2, y2 = p2
    
    f1 = ls(x1, y1)
    f2 = ls(x2, y2)
    
    # Check that f1 and f2 have opposite signs
    if f1 * f2 > 0
        error("Points don't bracket a zero of the function")
    end
    
    # Handle edge cases
    if abs(f1) < 1e-10
        return p1
    elseif abs(f2) < 1e-10
        return p2
    elseif abs(f1 - f2) < 1e-10
        return (0.5 * (x1 + x2), 0.5 * (y1 + y2))
    end
    
    # Linear interpolation parameter
    t = f1 / (f1 - f2)
    
    # Interpolated point
    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    
    return (x, y)
end

"""
    polygon_area(vertices)

Calcule l'aire d'un polygone défini par ses sommets.
"""
function polygon_area(vertices)
    n = length(vertices)
    if n < 3
        return 0.0
    end
    
    area = 0.0
    for i in 1:n
        j = i % n + 1
        area += vertices[i][1] * vertices[j][2]
        area -= vertices[j][1] * vertices[i][2]
    end
    
    return 0.5 * abs(area)
end

# Add this helper function to calculate the centroid of a polygon
function polygon_centroid(vertices)
    n = length(vertices)
    if n < 3
        # Not a proper polygon, return center
        x_sum = sum(v[1] for v in vertices)
        y_sum = sum(v[2] for v in vertices)
        return (x_sum / n, y_sum / n)
    end
    
    area = 0.0
    cx = 0.0
    cy = 0.0
    
    for i in 1:n
        j = i % n + 1
        xi, yi = vertices[i]
        xj, yj = vertices[j]
        
        # Cross product term
        cross_term = xi*yj - xj*yi
        
        # Area contribution
        area += cross_term
        
        # Centroid contribution
        cx += (xi + xj) * cross_term
        cy += (yi + yj) * cross_term
    end
    
    # Complete the calculations
    area = 0.5 * abs(area)
    
    if area < 1e-10
        # Avoid division by near-zero
        x_sum = sum(v[1] for v in vertices)
        y_sum = sum(v[2] for v in vertices)
        return (x_sum / n, y_sum / n)
    end
    
    # The 1/6A factor: 1/6 from the formula and 1/2 from the area calculation
    factor = 1.0 / (6.0 * area)
    cx = cx * factor
    cy = cy * factor
    
    return (abs(cx), abs(cy))
end

"""
    getcc(ls, xcell, ycell; nsub=0)

Calcule la fraction volumique (en 2D: fraction d'aire) de la cellule occupée par le fluide (ls < 0).
"""
function getcc(ls, xcell, ycell; nsub=0)
    # Test rapide du type de cellule
    ctype = getcelltype(ls, xcell, ycell)
    if ctype == 1
        return 1.0  # cellule pleine
    elseif ctype == 0
        return 0.0  # cellule vide
    end
    
    x_min, x_max = xcell
    y_min, y_max = ycell
    dx = x_max - x_min
    dy = y_max - y_min
    
    # Les 4 sommets de la cellule
    corners = [
        (x_min, y_min),  # coin bas-gauche
        (x_max, y_min),  # coin bas-droite
        (x_max, y_max),  # coin haut-droite
        (x_min, y_max)   # coin haut-gauche
    ]
    
    # Évaluer la level set aux 4 sommets
    corner_values = [ls(x, y) for (x, y) in corners]
    
    # Sommets intérieurs (ls < 0)
    interior_corners = [corners[i] for i in 1:4 if corner_values[i] < 0]
    
    # Arêtes de la cellule
    edges = [
        (corners[1], corners[2]),  # bas
        (corners[2], corners[3]),  # droite
        (corners[3], corners[4]),  # haut
        (corners[4], corners[1])   # gauche
    ]
    
    # Trouver les intersections sur les arêtes
    intersections = []
    for (p1, p2) in edges
        f1 = ls(p1[1], p1[2])
        f2 = ls(p2[1], p2[2])
        
        if f1 * f2 <= 0 && !(abs(f1) < 1e-10 && abs(f2) < 1e-10)
            # L'arête traverse l'interface
            if abs(f1) < 1e-10
                push!(intersections, p1)
            elseif abs(f2) < 1e-10
                push!(intersections, p2)
            else
                # Trouver le point d'intersection précisément
                intersection = find_zero(ls, p1, p2)
                push!(intersections, intersection)
            end
        end
    end
    
    # Construire le polygone (intérieur + intersections)
    if isempty(interior_corners)
        # Cas spécial: aucun coin intérieur, mais cellule coupée
        # On utilise juste les intersections comme polygone
        if length(intersections) >= 3
            # Trier les intersections pour former un polygone valide
            center_x = sum(p[1] for p in intersections) / length(intersections)
            center_y = sum(p[2] for p in intersections) / length(intersections)
            
            sort!(intersections, by = p -> atan(p[2] - center_y, p[1] - center_x))
            
            # Calculer l'aire du polygone et normaliser
            wet_area = polygon_area(intersections)
            return wet_area / (dx * dy)
        else
            # Moins de 3 intersections, utiliser approximation
            return 0.0
        end
    else
        # Cas général: combiner les coins intérieurs et les intersections
        vertices = vcat(interior_corners, intersections)
        
        # Trier les sommets pour former un polygone valide
        center_x = sum(p[1] for p in vertices) / length(vertices)
        center_y = sum(p[2] for p in vertices) / length(vertices)
        
        sort!(vertices, by = p -> atan(p[2] - center_y, p[1] - center_x))
        
        # Calculer l'aire du polygone et normaliser
        wet_area = polygon_area(vertices)
        return wet_area / (dx * dy)
    end
end

"""
    compute_interface_length(ls, xcell, ycell)

Compute the length of the interface segment within a 2D cell.
Returns the interface length for cut cells, or 0.0 for full/empty cells.
"""
function compute_interface_length(ls, xcell, ycell)
    # First check if there's an interface in the cell
    ctype = getcelltype(ls, xcell, ycell)
    if ctype != -1  # Not a cut cell
        return 0.0
    end
    
    x_min, x_max = xcell
    y_min, y_max = ycell
    
    # The 4 corners of the cell
    corners = [
        (x_min, y_min),  # bottom-left
        (x_max, y_min),  # bottom-right
        (x_max, y_max),  # top-right
        (x_min, y_max)   # top-left
    ]
    
    # The 4 edges of the cell
    edges = [
        (corners[1], corners[2]),  # bottom
        (corners[2], corners[3]),  # right
        (corners[3], corners[4]),  # top
        (corners[4], corners[1])   # left
    ]
    
    # Find intersections on edges
    intersections = []
    for (p1, p2) in edges
        f1 = ls(p1[1], p1[2])
        f2 = ls(p2[1], p2[2])
        
        if f1 * f2 <= 0 && !(abs(f1) < 1e-10 && abs(f2) < 1e-10)
            # Edge crosses interface
            if abs(f1) < 1e-10
                push!(intersections, p1)
            elseif abs(f2) < 1e-10
                push!(intersections, p2)
            else
                # Find intersection point using linear interpolation
                intersection = find_zero(ls, p1, p2)
                push!(intersections, intersection)
            end
        end
    end
    
    # Standard case: exactly 2 intersection points
    if length(intersections) == 2
        p1, p2 = intersections
        # Compute Euclidean distance between the two intersection points
        return sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    elseif length(intersections) > 2
        # Complex case with multiple intersections
        # We need to order the intersections to calculate the total interface length
        
        # Group intersections by edge
        edge_intersections = [[] for _ in 1:4]
        for p in intersections
            for (i, (p1, p2)) in enumerate(edges)
                # Check if point is on this edge using distance calculations
                on_edge = is_point_on_segment(p, p1, p2)
                if on_edge
                    push!(edge_intersections[i], p)
                    break
                end
            end
        end
        
        # Get ordered list of intersections following the cell perimeter
        ordered_intersections = []
        for i in 1:4
            if !isempty(edge_intersections[i])
                # Sort points along this edge
                p1, p2 = edges[i]
                if abs(p2[1] - p1[1]) > abs(p2[2] - p1[2])
                    # Edge is more horizontal, sort by x
                    sort!(edge_intersections[i], by = p -> p[1])
                else
                    # Edge is more vertical, sort by y
                    sort!(edge_intersections[i], by = p -> p[2])
                end
                append!(ordered_intersections, edge_intersections[i])
            end
        end
        
        # Calculate total interface length by connecting alternate pairs
        interface_length = 0.0
        for i in 1:2:length(ordered_intersections)-1
            p1 = ordered_intersections[i]
            p2 = ordered_intersections[i+1]
            interface_length += sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        end
        
        return interface_length
    else
        # Less than 2 intersections (shouldn't happen for valid cut cell)
        return 0.0
    end
end

# Helper function to check if a point is on a line segment
function is_point_on_segment(p, p1, p2, tol=1e-10)
    # Calculate distances
    d_total = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    d1 = sqrt((p[1] - p1[1])^2 + (p[2] - p1[2])^2)
    d2 = sqrt((p[1] - p2[1])^2 + (p[2] - p2[2])^2)
    
    # Point is on segment if sum of distances equals total length (with tolerance)
    return abs(d1 + d2 - d_total) < tol
end

# 3D functions
"""
    getcelltype(ls, xcell, ycell, zcell)
    
Returns cell type for a 3D cell:
- 1: full cell (level set < 0 everywhere)
- 0: empty cell (level set > 0 everywhere)
- -1: cut cell (level set changes sign)
"""
function getcelltype(ls, xcell, ycell, zcell)
    # Evaluate level set at the 8 vertices of the cell
    φ = [
        ls(xcell[1], ycell[1], zcell[1]),  # bottom-left-back
        ls(xcell[2], ycell[1], zcell[1]),  # bottom-right-back
        ls(xcell[2], ycell[2], zcell[1]),  # top-right-back
        ls(xcell[1], ycell[2], zcell[1]),  # top-left-back
        ls(xcell[1], ycell[1], zcell[2]),  # bottom-left-front
        ls(xcell[2], ycell[1], zcell[2]),  # bottom-right-front
        ls(xcell[2], ycell[2], zcell[2]),  # top-right-front
        ls(xcell[1], ycell[2], zcell[2])   # top-left-front
    ]
    
    if all(φ .< 0)
        return 1  # full cell
    elseif all(φ .> 0)
        return 0  # empty cell
    else
        return -1  # cut cell
    end
end

"""
    find_zero_3d(ls, p1, p2)
    
Finds the zero of the level set function on a 3D segment [p1, p2] using linear interpolation.
"""
function find_zero_3d(ls, p1, p2)
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    
    f1 = ls(x1, y1, z1)
    f2 = ls(x2, y2, z2)
    
    # Check that f1 and f2 have opposite signs
    if f1 * f2 > 0
        error("Points don't bracket a zero of the function")
    end
    
    # Handle edge cases
    if abs(f1) < 1e-10
        return p1
    elseif abs(f2) < 1e-10
        return p2
    elseif abs(f1 - f2) < 1e-10
        return (0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (z1 + z2))
    end
    
    # Linear interpolation parameter
    t = f1 / (f1 - f2)
    
    # Interpolated point
    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    z = z1 + t * (z2 - z1)
    
    return (x, y, z)
end

"""
    tetrahedron_volume(v0, v1, v2, v3)
    
Compute the signed volume of a tetrahedron defined by four vertices.
"""
function tetrahedron_volume(v0, v1, v2, v3)
    # Create vectors from v0 to other vertices
    a = (v1[1] - v0[1], v1[2] - v0[2], v1[3] - v0[3])
    b = (v2[1] - v0[1], v2[2] - v0[2], v2[3] - v0[3])
    c = (v3[1] - v0[1], v3[2] - v0[2], v3[3] - v0[3])
    
    # Triple scalar product (v1-v0) ⋅ ((v2-v0) × (v3-v0))
    volume = (
        a[1] * (b[2] * c[3] - b[3] * c[2]) +
        a[2] * (b[3] * c[1] - b[1] * c[3]) +
        a[3] * (b[1] * c[2] - b[2] * c[1])
    ) / 6.0
    
    return abs(volume)  # We take absolute value as we're concerned with magnitude
end

"""
    getcc(ls, xcell, ycell, zcell; nsub=0)
    
Computes the volume fraction of a 3D cell occupied by fluid (ls < 0)
using the generalized shoelace formula with an accurate triangulation approach.
"""
function getcc(ls, xcell, ycell, zcell; nsub=0)
    # Check trivial cases
    ctype = getcelltype(ls, xcell, ycell, zcell)
    if ctype == 1
        return 1.0  # full cell
    elseif ctype == 0
        return 0.0  # empty cell
    end
    
    # Cell dimensions
    x_min, x_max = xcell
    y_min, y_max = ycell
    z_min, z_max = zcell
    
    cell_volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
    
    # The 8 vertices of the cell
    corners = [
        (x_min, y_min, z_min),  # 1: bottom-left-back
        (x_max, y_min, z_min),  # 2: bottom-right-back
        (x_max, y_max, z_min),  # 3: top-right-back
        (x_min, y_max, z_min),  # 4: top-left-back
        (x_min, y_min, z_max),  # 5: bottom-left-front
        (x_max, y_min, z_max),  # 6: bottom-right-front
        (x_max, y_max, z_max),  # 7: top-right-front
        (x_min, y_max, z_max)   # 8: top-left-front
    ]
    
    # Evaluate level set at corners
    corner_values = [ls(x, y, z) for (x, y, z) in corners]
    
    # Find interior corners (ls < 0)
    interior_indices = findall(v -> v < 0, corner_values)
    interior_corners = [corners[i] for i in interior_indices]
    
    # The 12 edges of the cell
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 1),  # bottom face edges
        (5, 6), (6, 7), (7, 8), (8, 5),  # top face edges
        (1, 5), (2, 6), (3, 7), (4, 8)   # vertical edges
    ]
    
    # Find intersections on edges
    intersection_points = []
    for (i, j) in edges
        p1, p2 = corners[i], corners[j]
        f1, f2 = corner_values[i], corner_values[j]
        
        if f1 * f2 < 0
            # Edge crosses the interface, find intersection
            intersection = find_zero_3d(ls, p1, p2)
            push!(intersection_points, intersection)
        elseif abs(f1) < 1e-10
            # Vertex is on interface
            push!(intersection_points, p1)
        elseif abs(f2) < 1e-10
            # Vertex is on interface
            push!(intersection_points, p2)
        end
    end
    
    # Remove duplicate intersection points
    unique_intersections = []
    for p in intersection_points
        is_duplicate = false
        for u in unique_intersections
            if sqrt((p[1] - u[1])^2 + (p[2] - u[2])^2 + (p[3] - u[3])^2) < 1e-10
                is_duplicate = true
                break
            end
        end
        if !is_duplicate
            push!(unique_intersections, p)
        end
    end
    intersection_points = unique_intersections
    
    # Based on the number of interior corners, handle different cases
    interior_count = length(interior_indices)

    if interior_count == 0
        # No interior corners - the interface cuts the cell forming a small polyhedron
        # We need to properly triangulate the interface polygon
        if length(intersection_points) < 3
            return 0.0
        end
        
        # Find centroid of intersection points to use as reference
        centroid = (
            sum(p[1] for p in intersection_points) / length(intersection_points),
            sum(p[2] for p in intersection_points) / length(intersection_points),
            sum(p[3] for p in intersection_points) / length(intersection_points)
        )
        
        # Orient the points to establish a consistent ordering
        # Project each point onto a unit sphere centered at the centroid
        # Then sort based on spherical coordinates (theta, phi)
        function get_spherical_coords(p)
            dx = p[1] - centroid[1]
            dy = p[2] - centroid[2]
            dz = p[3] - centroid[3]
            r = sqrt(dx^2 + dy^2 + dz^2)
            if r < 1e-10
                return (0.0, 0.0)
            end
            theta = atan(dy, dx)
            phi = acos(dz / r)
            return (theta, phi)
        end
        
        # First, detect if the points are roughly planar
        if length(intersection_points) >= 4
            # Check if the points lie approximately in a plane
            # by comparing the determinants of the point triplets
            is_planar = true
            p1, p2, p3 = intersection_points[1:3]
            
            # Create vectors from p1 to p2 and p1 to p3
            v1 = (p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
            v2 = (p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])
            
            # Cross product to find the normal
            normal = (
                v1[2]*v2[3] - v1[3]*v2[2],
                v1[3]*v2[1] - v1[1]*v2[3],
                v1[1]*v2[2] - v1[2]*v2[1]
            )
            
            # Normalize the normal vector
            norm_length = sqrt(normal[1]^2 + normal[2]^2 + normal[3]^2)
            if norm_length > 1e-10
                normal = (normal[1]/norm_length, normal[2]/norm_length, normal[3]/norm_length)
            end
            
            # Use the centroid to orient the normal toward the fluid region
            if ls(centroid...) > 0
                normal = (-normal[1], -normal[2], -normal[3])
            end
            
            # For planar interfaces, we can use a 2D triangulation approach
            # Project points onto the plane and use a 2D ordering
            function project_to_2d(p)
                # Find the dominant axis of the normal
                _, dominant_axis = findmax([abs(normal[1]), abs(normal[2]), abs(normal[3])])
                
                if dominant_axis == 1
                    # X is dominant, project to YZ plane
                    return (p[2], p[3])
                elseif dominant_axis == 2
                    # Y is dominant, project to XZ plane
                    return (p[1], p[3])
                else
                    # Z is dominant, project to XY plane
                    return (p[1], p[2])
                end
            end
            
            # Project points to 2D
            projected_points = [project_to_2d(p) for p in intersection_points]
            
            # Find 2D centroid
            centroid_2d = (
                sum(p[1] for p in projected_points) / length(projected_points),
                sum(p[2] for p in projected_points) / length(projected_points)
            )
            
            # Sort points based on angle in 2D plane
            sorted_indices = sortperm(projected_points, by = p -> 
                atan(p[2] - centroid_2d[2], p[1] - centroid_2d[1]))
            sorted_points = intersection_points[sorted_indices]
        else
            # For non-planar or small number of points, use 3D approach
            sorted_points = sort(intersection_points, by = get_spherical_coords)
        end
        
        # Create triangulation using the sorted points
        # Each triangle connects the centroid to a pair of consecutive points
        triangles = []
        for i in 1:length(sorted_points)
            j = (i % length(sorted_points)) + 1
            push!(triangles, (centroid, sorted_points[i], sorted_points[j]))
        end
        
        # Calculate volume using generalized shoelace formula
        volume = 0.0
        for (v1, v2, v3) in triangles
            # Create a tetrahedron with the origin and compute its volume
            origin = (0.0, 0.0, 0.0)
            volume += tetrahedron_volume(origin, v1, v2, v3)
        end
        
        return volume / cell_volume
        
    elseif interior_count == 8
        # All corners are inside
        return 1.0
        
    else
        # Some corners are inside, some outside
        # We'll use the interface to split the cell into inside and outside parts
        
        # Sort the intersection points to create a valid triangulation of the interface
        if length(intersection_points) >= 3
            # Find centroid of the interface points
            interface_centroid = (
                sum(p[1] for p in intersection_points) / length(intersection_points),
                sum(p[2] for p in intersection_points) / length(intersection_points),
                sum(p[3] for p in intersection_points) / length(intersection_points)
            )
            
            # Compute a normal vector for the interface points
            plane_points = intersection_points[1:min(3, length(intersection_points))]
            v1 = (plane_points[2][1] - plane_points[1][1],
                plane_points[2][2] - plane_points[1][2],
                plane_points[2][3] - plane_points[1][3])
            
            v2 = (plane_points[3 % length(plane_points) + 1][1] - plane_points[1][1],
                plane_points[3 % length(plane_points) + 1][2] - plane_points[1][2],
                plane_points[3 % length(plane_points) + 1][3] - plane_points[1][3])
            
            normal = (
                v1[2]*v2[3] - v1[3]*v2[2],
                v1[3]*v2[1] - v1[1]*v2[3],
                v1[1]*v2[2] - v1[2]*v2[1]
            )
            
            # Normalize
            norm_length = sqrt(normal[1]^2 + normal[2]^2 + normal[3]^2)
            if norm_length > 1e-10
                normal = (normal[1]/norm_length, normal[2]/norm_length, normal[3]/norm_length)
            else
                # Fallback if we can't compute a good normal
                normal = (0.0, 0.0, 1.0)
            end
            
            # Create a coordinate system on the plane
            u_axis = (0.0, 0.0, 0.0)
            if abs(normal[3]) > 0.9
                # If normal is close to z-axis, use x-axis for u
                u_axis = (1.0, 0.0, 0.0)
            else
                # Otherwise use z-axis
                u_axis = (0.0, 0.0, 1.0)
            end
            
            # v = normal × u
            v_axis = (
                normal[2]*u_axis[3] - normal[3]*u_axis[2],
                normal[3]*u_axis[1] - normal[1]*u_axis[3],
                normal[1]*u_axis[2] - normal[2]*u_axis[1]
            )
            
            # Ensure orthogonality: u = v × normal
            u_axis = (
                v_axis[2]*normal[3] - v_axis[3]*normal[2],
                v_axis[3]*normal[1] - v_axis[1]*normal[3],
                v_axis[1]*normal[2] - v_axis[2]*normal[1]
            )
            
            # Normalize u and v
            u_length = sqrt(u_axis[1]^2 + u_axis[2]^2 + u_axis[3]^2)
            u_axis = (u_axis[1]/u_length, u_axis[2]/u_length, u_axis[3]/u_length)
            
            v_length = sqrt(v_axis[1]^2 + v_axis[2]^2 + v_axis[3]^2)
            v_axis = (v_axis[1]/v_length, v_axis[2]/v_length, v_axis[3]/v_length)
            
            # Project points to get 2D coordinates
            function project_to_plane(p)
                # Vector from centroid to point
                dp = (p[1] - interface_centroid[1],
                    p[2] - interface_centroid[2],
                    p[3] - interface_centroid[3])
                
                # Projections onto u and v axes
                u_coord = dp[1]*u_axis[1] + dp[2]*u_axis[2] + dp[3]*u_axis[3]
                v_coord = dp[1]*v_axis[1] + dp[2]*v_axis[2] + dp[3]*v_axis[3]
                
                return (u_coord, v_coord)
            end
            
            # Project and sort
            projected_points = [project_to_plane(p) for p in intersection_points]
            
            # Sort by angle from the origin in the projected 2D plane
            sorted_indices = sortperm(projected_points, by = p -> atan(p[2], p[1]))
            sorted_points = intersection_points[sorted_indices]
            
            # Create triangles for the interface
            interface_triangles = []
            for i in 2:(length(sorted_points)-1)
                push!(interface_triangles, (sorted_points[1], sorted_points[i], sorted_points[i+1]))
            end
            
            # Calculate volume by summing tetrahedra formed between each interface triangle
            # and each interior corner
            volume = 0.0
            
            # Form tetrahedra between interior corners and interface triangles
            for corner in interior_corners
                for (v1, v2, v3) in interface_triangles
                    volume += tetrahedron_volume(corner, v1, v2, v3)
                end
            end
            
            # Add volumes for any cell faces that are completely inside the fluid
            # This handles cases where the interface doesn't cut through all faces
            cell_faces = [
                [(1,2,3), (1,3,4)],   # bottom face (-z)
                [(5,6,7), (5,7,8)],   # top face (+z)
                [(1,2,6), (1,6,5)],   # front face (-y)
                [(3,4,8), (3,8,7)],   # back face (+y)
                [(1,4,8), (1,8,5)],   # left face (-x)
                [(2,3,7), (2,7,6)]    # right face (+x)
            ]
            
            # Reference point for each face (used to check if face is inside fluid)
            face_centers = [
                ((x_min + x_max)/2, (y_min + y_max)/2, z_min),           # bottom
                ((x_min + x_max)/2, (y_min + y_max)/2, z_max),           # top
                ((x_min + x_max)/2, y_min, (z_min + z_max)/2),           # front
                ((x_min + x_max)/2, y_max, (z_min + z_max)/2),           # back
                (x_min, (y_min + y_max)/2, (z_min + z_max)/2),           # left
                (x_max, (y_min + y_max)/2, (z_min + z_max)/2)            # right
            ]
            
            # Check each face of the cell
            for (face_idx, face) in enumerate(cell_faces)
                # If the face center is inside the fluid, add the face's contribution
                if ls(face_centers[face_idx]...) < 0
                    is_intact = true
                    for (i, j, k) in face
                        # If any vertex of this face triangle is outside, the face is cut
                        if any(v -> v > 0, [corner_values[i], corner_values[j], corner_values[k]])
                            is_intact = false
                            break
                        end
                    end
                    
                    # If the face is intact inside the fluid, add a pyramid from centroid
                    if is_intact
                        # Add the volume contribution from this face
                        for (i, j, k) in face
                            v1, v2, v3 = corners[i], corners[j], corners[k]
                            face_centroid = (
                                (v1[1] + v2[1] + v3[1])/3,
                                (v1[2] + v2[2] + v3[2])/3,
                                (v1[3] + v2[3] + v3[3])/3
                            )
                            # Form tetrahedron with cell center
                            cell_center = (
                                (x_min + x_max)/2,
                                (y_min + y_max)/2,
                                (z_min + z_max)/2
                            )
                            volume += tetrahedron_volume(cell_center, v1, v2, v3)
                        end
                    end
                end
            end
            
            # Return normalized volume fraction, ensuring it's in [0,1]
            return min(1.0, max(0.0, volume / cell_volume))
        else
            # Not enough intersection points for a valid interface
            # Use the proportion of interior corners as an approximation
            return interior_count / 8.0
        end
    end
end

"""
    compute_interface_area(ls, xcell, ycell, zcell)

Compute the area of the interface within a 3D cell.
Returns the interface area for cut cells, or 0.0 for full/empty cells.
"""
function compute_interface_area(ls, xcell, ycell, zcell)
    # First check if there's an interface in the cell
    ctype = getcelltype(ls, xcell, ycell, zcell)
    if ctype != -1  # Not a cut cell
        return 0.0
    end
    
    # Cell dimensions
    x_min, x_max = xcell
    y_min, y_max = ycell
    z_min, z_max = zcell
    
    # The 8 vertices of the cell
    corners = [
        (x_min, y_min, z_min),  # 1: bottom-left-back
        (x_max, y_min, z_min),  # 2: bottom-right-back
        (x_max, y_max, z_min),  # 3: top-right-back
        (x_min, y_max, z_min),  # 4: top-left-back
        (x_min, y_min, z_max),  # 5: bottom-left-front
        (x_max, y_min, z_max),  # 6: bottom-right-front
        (x_max, y_max, z_max),  # 7: top-right-front
        (x_min, y_max, z_max)   # 8: top-left-front
    ]
    
    # Evaluate level set at corners
    corner_values = [ls(x, y, z) for (x, y, z) in corners]
    
    # The 12 edges of the cell
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 1),  # bottom face edges
        (5, 6), (6, 7), (7, 8), (8, 5),  # top face edges
        (1, 5), (2, 6), (3, 7), (4, 8)   # vertical edges
    ]
    
    # Find intersections on edges
    intersection_points = []
    for (i, j) in edges
        p1, p2 = corners[i], corners[j]
        f1, f2 = corner_values[i], corner_values[j]
        
        if f1 * f2 < 0
            # Edge crosses the interface, find intersection
            intersection = find_zero_3d(ls, p1, p2)
            push!(intersection_points, intersection)
        elseif abs(f1) < 1e-10
            # Vertex is on interface
            push!(intersection_points, p1)
        elseif abs(f2) < 1e-10
            # Vertex is on interface
            push!(intersection_points, p2)
        end
    end
    
    # Remove duplicate intersection points
    unique_intersections = []
    for p in intersection_points
        is_duplicate = false
        for u in unique_intersections
            if sqrt((p[1] - u[1])^2 + (p[2] - u[2])^2 + (p[3] - u[3])^2) < 1e-10
                is_duplicate = true
                break
            end
        end
        if !is_duplicate
            push!(unique_intersections, p)
        end
    end
    intersection_points = unique_intersections
    
    # Not enough points to form a valid interface
    if length(intersection_points) < 3
        return 0.0
    end
    
    # Find centroid of intersection points
    centroid = (
        sum(p[1] for p in intersection_points) / length(intersection_points),
        sum(p[2] for p in intersection_points) / length(intersection_points),
        sum(p[3] for p in intersection_points) / length(intersection_points)
    )
    
    # Compute a normal vector for the interface points
    plane_points = intersection_points[1:min(3, length(intersection_points))]
    v1 = (plane_points[2][1] - plane_points[1][1],
          plane_points[2][2] - plane_points[1][2],
          plane_points[2][3] - plane_points[1][3])
    
    v2 = (plane_points[3 % length(plane_points) + 1][1] - plane_points[1][1],
          plane_points[3 % length(plane_points) + 1][2] - plane_points[1][2],
          plane_points[3 % length(plane_points) + 1][3] - plane_points[1][3])
    
    normal = (
        v1[2]*v2[3] - v1[3]*v2[2],
        v1[3]*v2[1] - v1[1]*v2[3],
        v1[1]*v2[2] - v1[2]*v2[1]
    )
    
    # Normalize
    norm_length = sqrt(normal[1]^2 + normal[2]^2 + normal[3]^2)
    if norm_length > 1e-10
        normal = (normal[1]/norm_length, normal[2]/norm_length, normal[3]/norm_length)
    else
        # Fallback if we can't compute a good normal
        normal = (0.0, 0.0, 1.0)
    end
    
    # Check if the points are roughly planar
    if length(intersection_points) >= 4
        # For planar interfaces, create a 2D coordinate system on the plane
        u_axis = (0.0, 0.0, 0.0)
        if abs(normal[3]) > 0.9
            # If normal is close to z-axis, use x-axis for u
            u_axis = (1.0, 0.0, 0.0)
        else
            # Otherwise use z-axis
            u_axis = (0.0, 0.0, 1.0)
        end
        
        # v = normal × u
        v_axis = (
            normal[2]*u_axis[3] - normal[3]*u_axis[2],
            normal[3]*u_axis[1] - normal[1]*u_axis[3],
            normal[1]*u_axis[2] - normal[2]*u_axis[1]
        )
        
        # Now recompute u to ensure orthogonality: u = v × normal
        u_axis = (
            v_axis[2]*normal[3] - v_axis[3]*normal[2],
            v_axis[3]*normal[1] - v_axis[1]*normal[3],
            v_axis[1]*normal[2] - v_axis[2]*normal[1]
        )
        
        # Normalize u and v
        u_length = sqrt(u_axis[1]^2 + u_axis[2]^2 + u_axis[3]^2)
        u_axis = (u_axis[1]/u_length, u_axis[2]/u_length, u_axis[3]/u_length)
        
        v_length = sqrt(v_axis[1]^2 + v_axis[2]^2 + v_axis[3]^2)
        v_axis = (v_axis[1]/v_length, v_axis[2]/v_length, v_axis[3]/v_length)
        
        # Project points to get 2D coordinates
        function project_to_plane(p)
            # Vector from centroid to point
            dp = (p[1] - centroid[1],
                  p[2] - centroid[2],
                  p[3] - centroid[3])
            
            # Projections onto u and v axes
            u_coord = dp[1]*u_axis[1] + dp[2]*u_axis[2] + dp[3]*u_axis[3]
            v_coord = dp[1]*v_axis[1] + dp[2]*v_axis[2] + dp[3]*v_axis[3]
            
            return (u_coord, v_coord)
        end
        
        # Project and sort
        projected_points = [project_to_plane(p) for p in intersection_points]
        
        # Sort by angle from the origin in the projected 2D plane
        sorted_indices = sortperm(projected_points, by = p -> atan(p[2], p[1]))
        sorted_points = intersection_points[sorted_indices]
        
        # Create triangulation using fan triangulation from the first point
        triangles = []
        for i in 2:(length(sorted_points)-1)
            push!(triangles, (sorted_points[1], sorted_points[i], sorted_points[i+1]))
        end
        
        # Calculate the total area of the interface
        area = 0.0
        for (v1, v2, v3) in triangles
            # Calculate area of this triangle using the cross product formula
            e1 = (v2[1] - v1[1], v2[2] - v1[2], v2[3] - v1[3])
            e2 = (v3[1] - v1[1], v3[2] - v1[2], v3[3] - v1[3])
            
            # Cross product
            cross = (
                e1[2]*e2[3] - e1[3]*e2[2],
                e1[3]*e2[1] - e1[1]*e2[3],
                e1[1]*e2[2] - e1[2]*e2[1]
            )
            
            # Area = 1/2 * |cross|
            area += 0.5 * sqrt(cross[1]^2 + cross[2]^2 + cross[3]^2)
        end
        
        return area
    else
        # Not enough points for a good triangulation
        return 0.0
    end
end


# Extensions pour 4D
# TODO: Implémenter les méthodes getcelltype et getcc pour  4D


end # module