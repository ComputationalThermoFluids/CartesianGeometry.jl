function vofinit!(xex, f, x::Number)
    t = f(x)

    val = one(x)

    isnonpositive(t) && return val
    isnonnegative(t) && return zero(val)

    nothing
end

"""

1d implementation complies with Vofi 2.0's API (see 2d and 3d).

!!! note "Convention"

    Even if cell is empty, return full centroid coordinates.

"""
function vofinit!(xex, f, x::SVector; nex=Cint.((1, 1)))
    t = SVector{2}(f(i) for i in x)

    val = x[2] - x[1]

    if all(isnonpositive, t)
        isone(first(nex)) && (xex[1] = sum(x) / 2)
        isone(last(nex)) && (xex[end] = zero(xex[1]))
        return val
    end

    if all(isnonnegative, t)
        isone(first(nex)) && (xex[1] = sum(x) / 2)
        isone(last(nex)) && (xex[end] = zero(xex[1]))
        return zero(val)
    end

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    if isnonnegative(t[1])
        isone(first(nex)) && (xex[1] = (ξ + x[2]) / 2)
        isone(last(nex)) && (xex[end] = 1.0)
        return x[2] - ξ
    end

    if isnonnegative(t[2])
        isone(first(nex)) && (xex[1] = (x[1] + ξ) / 2)
        isone(last(nex)) && (xex[end] = 1.0)
        return ξ - x[1]
    end

    nothing
end

function vofinit!(xex, f, x::SVector{2}, y::Number)
    t = SVector{2}(f(i, y) for i in x)

    val = x[2] - x[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return x[2] - ξ
    isnonnegative(t[2]) && return ξ - x[1]

    nothing
end

function vofinit!(xex, f, x::Number, y::SVector{2})
    t = SVector{2}(f(x, j) for j in y)

    val = y[2] - y[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (y[2] * t[1] - y[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return y[2] - ξ
    isnonnegative(t[2]) && return ξ - y[1]

    nothing
end

"""

Call Vofi 2.0 for exact integration. ###

"""
function vofinit!(xex, f, x::SVector, y::SVector; nex=Cint.((1, 1)))
    t = SMatrix{2,2}(f(i, j) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[end] = zero(xex[1])
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[end] = zero(xex[1])
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    val * getcc(f, x0, h0, xex, Cint(2); nex)
end

function vofinit!(xex, f, x::Number, y::SVector{2}, z::SVector{2})
    t = SMatrix{2,2}(f(x, j, k) for j in y, k in z)

    val = (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((y[1], z[1], zero(val)))
    h0 = Cdouble.((y[2]-y[1], z[2]-z[1], one(val)))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(2); nex) do y, z, _
        f(x, y, z)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::Number, z::SVector{2})
    t = SMatrix{2,2}(f(i, y, k) for i in x, k in z)

    val = (z[2] - z[1]) * (x[2] - x[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((z[1], x[1], zero(val)))
    h0 = Cdouble.((z[2]-z[1], x[2]-x[1], one(val)))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(2); nex) do z, x, _
        f(x, y, z)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::Number)
    t = SMatrix{2,2}(f(i, j, z) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(2); nex) do x, y, _
        f(x, y, z)
    end
end

"""

Call Vofi 2.0 for exact integration.

"""
function vofinit!(xex, f, x::SVector, y::SVector, z::SVector; nex=Cint.((1, 1)))
    t = SArray{Tuple{2,2,2}}(f(i, j, k) for i in x, j in y, k in z)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[end] = zero(xex[1])
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[end] = zero(xex[1])
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    val * getcc(f, x0, h0, xex, Cint(3); nex)
end

"""

Call Vofi extended for exact 4D integration.

"""
function vofinit!(xex, f, x::SVector, y::SVector, z::SVector, w::SVector; nex=Cint.((1, 1)))
    t = SArray{Tuple{2,2,2,2}}(f(i, j, k, l) for i in x, j in y, k in z, l in w)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1]) * (w[2] - w[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[4] = sum(w) / 2
        xex[end] = zero(xex[1])
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[4] = sum(w) / 2
        xex[end] = zero(xex[1])
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], z[1], w[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1], w[2]-w[1]))

    val * getcc(f, x0, h0, xex, Cint(4); nex)
end

# Fonctions pour les combinaisons d'arguments 4D avec une dimension fixe
function vofinit!(xex, f, x::Number, y::SVector{2}, z::SVector{2}, w::SVector{2})
    t = SArray{Tuple{2,2,2}}(f(x, j, k, l) for j in y, k in z, l in w)

    val = (y[2] - y[1]) * (z[2] - z[1]) * (w[2] - w[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((y[1], z[1], w[1]))
    h0 = Cdouble.((y[2]-y[1], z[2]-z[1], w[2]-w[1]))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(3); nex) do y, z, w
        f(x, y, z, w)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::Number, z::SVector{2}, w::SVector{2})
    t = SArray{Tuple{2,2,2}}(f(i, y, k, l) for i in x, k in z, l in w)

    val = (x[2] - x[1]) * (z[2] - z[1]) * (w[2] - w[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], z[1], w[1]))
    h0 = Cdouble.((x[2]-x[1], z[2]-z[1], w[2]-w[1]))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(3); nex) do x, z, w
        f(x, y, z, w)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::Number, w::SVector{2})
    t = SArray{Tuple{2,2,2}}(f(i, j, z, l) for i in x, j in y, l in w)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (w[2] - w[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], w[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], w[2]-w[1]))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(3); nex) do x, y, w
        f(x, y, z, w)
    end
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::SVector{2}, w::Number)
    t = SArray{Tuple{2,2,2}}(f(i, j, k, w) for i in x, j in y, k in z)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    nex = Cint.((0, 0))
    val * getcc(x0, h0, xex, Cint(3); nex) do x, y, z
        f(x, y, z, w)
    end
end

function get_cell_type(f, x::SVector{2})
    t = SVector{2}(f(i) for i in x)
    if all(isnonpositive, t)
        return 1.0
    end
    if all(isnonnegative, t)
        return 0.0
    end
    if isnonnegative(t[1])
        return -1.0
    end
    if isnonnegative(t[2])
        return -1.0
    end
    nothing
end
function get_cell_type(f, x::SVector{2}, y::SVector{2})
    x0 = Cdouble.((x[1], y[1], 0.0))
    h0 = Cdouble.((x[2] - x[1], y[2] - y[1], 1.0))
    return getcelltype(f, x0, h0, Cint(2))
end

function get_cell_type(f, x::SVector{2}, y::SVector{2}, z::SVector{2})
    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2] - x[1], y[2] - y[1], z[2] - z[1]))
    return getcelltype(f, x0, h0, Cint(3))
end

function get_cell_type(f, x::SVector{2}, y::SVector{2}, z::SVector{2}, w::SVector{2})
    x0 = Cdouble.((x[1], y[1], z[1], w[1]))
    h0 = Cdouble.((x[2] - x[1], y[2] - y[1], z[2] - z[1], w[2] - w[1]))
    return getcelltype(f, x0, h0, Cint(4))
end



# Ajouter une fonction de dispatch à la fin du fichier:
"""
    vofinit_dispatch!(method, xex, f, args...)

Dispatche vers la méthode appropriée selon le paramètre 'method'.
"""
function vofinit_dispatch!(method, xex, f, args...; nex=Cint.((1, 1)))
    if method === :simple
        # Adapter les arguments pour SimpleVOF
        return simplevof_wrapper(xex, f, args...; nex=nex)
    else
        # Méthode par défaut: VOFI
        return vofinit!(xex, f, args...)
    end
end

"""
    simplevof_wrapper(xex, f, args...)

Wrapper for using SimpleVOF instead of VOFI for volume and surface fraction calculations.
Handles all the same argument patterns as vofinit!.
"""
function simplevof_wrapper(xex, f, args...; nex=Cint.((1, 1)))
    # 1D cases : OK
    # Case 1: 1D point
    if length(args) == 1 && args[1] isa Number
        x = args[1]
        t = f(x)
        
        val = one(x)  # Use one(x) for type consistency like VOFI
        
        if isnonpositive(t)
            isone(first(nex)) && (xex[1] = x)
            isone(last(nex)) && (xex[end] = zero(val))
            return val
        elseif isnonnegative(t)
            isone(first(nex)) && (xex[1] = x)
            isone(last(nex)) && (xex[end] = zero(val))
            return zero(val)
        else
            # This case shouldn't happen with proper isnonpositive/isnonnegative
            # but including for robustness
            isone(first(nex)) && (xex[1] = x)
            isone(last(nex)) && (xex[end] = one(val))
            return 0.5 * val  # approximation
        end
    
    # Case 2: 1D interval
    elseif length(args) == 1 && args[1] isa SVector
        x = args[1]
        # Create t vector using same approach as VOFI
        t = SVector{2}(f(i) for i in x)
        val = x[2] - x[1]
        
        # Handle trivial cases using VOFI's approach
        if all(isnonpositive, t)
            isone(first(nex)) && (xex[1] = sum(x) / 2)
            isone(last(nex)) && (xex[end] = zero(xex[1]))
            return val
        end
        
        if all(isnonnegative, t)
            isone(first(nex)) && (xex[1] = sum(x) / 2)
            isone(last(nex)) && (xex[end] = zero(xex[1]))
            return zero(val)
        end
        
        # Calculate zero-crossing using VOFI's formulation
        # This is mathematically equivalent to SimpleVOF's formula but matches VOFI's code
        ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])
        
        # Handle interface cases using VOFI's conditions
        if isnonnegative(t[1])
            isone(first(nex)) && (xex[1] = (ξ + x[2]) / 2)
            isone(last(nex)) && (xex[end] = 1.0)
            return x[2] - ξ
        end
        
        if isnonnegative(t[2])
            isone(first(nex)) && (xex[1] = (x[1] + ξ) / 2)
            isone(last(nex)) && (xex[end] = 1.0)
            return ξ - x[1]
        end

    
    # 2D cases : OK
    # Case 3: x interval, y point (vertical line)
    elseif length(args) == 2 && args[1] isa SVector{2} && args[2] isa Number
        x, y = args[1], args[2]
        # Function to evaluate along the line
        fline = (xt) -> f(xt, y)
        # Evaluate at endpoints
        t1, t2 = fline(x[1]), fline(x[2])
        
        # Handle trivial cases
        if isnonpositive(t1) && isnonpositive(t2)
            #xex[1] = (x[1] + x[2])/2
            #xex[2] = y
            #xex[end] = 0.0
            return x[2] - x[1]
        elseif isnonnegative(t1) && isnonnegative(t2)
            #xex[1] = (x[1] + x[2])/2
            #xex[2] = y
            #xex[end] = 0.0
            return 0.0
        end
        
        # For interface case, find approximate location
        ξ = x[1] - t1 * (x[2] - x[1]) / (t2 - t1)
        if t1 <= 0
            #xex[1] = (x[1] + ξ)/2
            #xex[2] = y
            #xex[end] = 1.0
            return ξ - x[1]
        else
            #xex[1] = (ξ + x[2])/2
            #xex[2] = y
            #xex[end] = 1.0
            return x[2] - ξ
        end
    
    # Case 4: x point, y interval (horizontal line)
    elseif length(args) == 2 && args[1] isa Number && args[2] isa SVector{2}
        x, y = args[1], args[2]
        # Function to evaluate along the line
        fline = (yt) -> f(x, yt)
        # Evaluate at endpoints
        t1, t2 = fline(y[1]), fline(y[2])
        
        # Handle trivial cases
        if isnonpositive(t1) && isnonpositive(t2)
            #xex[1] = x
            #xex[2] = (y[1] + y[2])/2
            #xex[end] = 0.0
            return y[2] - y[1]
        elseif isnonnegative(t1) && isnonnegative(t2)
            #xex[1] = x
            #xex[2] = (y[1] + y[2])/2
            #xex[end] = 0.0
            return 0.0
        end
        
        # For interface case, find approximate location
        ξ = y[1] - t1 * (y[2] - y[1]) / (t2 - t1)
        if t1 <= 0
            #xex[1] = x
            #xex[2] = (y[1] + ξ)/2
            #xex[end] = 1.0
            return ξ - y[1]
        else
            #xex[1] = x
            #xex[2] = (ξ + y[2])/2
            #xex[end] = 1.0
            return y[2] - ξ
        end
    
    # Case 5: 2D cell
    elseif length(args) == 2 && args[1] isa SVector && args[2] isa SVector
        x, y = args[1], args[2]
        
        # First check trivial cases
        cell_type = SimpleVOF.getcelltype(f, (x[1], x[2]), (y[1], y[2]))
        
        if cell_type == 1  # Full cell
            xex[1] = (x[1] + x[2])/2
            xex[2] = (y[1] + y[2])/2
            xex[end] = 0.0
            return (x[2] - x[1]) * (y[2] - y[1])
        elseif cell_type == 0  # Empty cell
            xex[1] = (x[1] + x[2])/2
            xex[2] = (y[1] + y[2])/2
            xex[end] = 0.0
            return 0.0
        end
        
        # For cut cells, compute the vertices of the wet polygon
        corners = [
            (x[1], y[1]),  # bottom-left
            (x[2], y[1]),  # bottom-right
            (x[2], y[2]),  # top-right
            (x[1], y[2])   # top-left
        ]
        
        # Evaluate level set at corners
        corner_values = [f(xc, yc) for (xc, yc) in corners]
        
        # Find interior corners (where level set < 0)
        interior_corners = [corners[i] for i in 1:4 if corner_values[i] < 0]
        
        # Find intersections on edges
        edges = [
            (corners[1], corners[2]),  # bottom
            (corners[2], corners[3]),  # right
            (corners[3], corners[4]),  # top
            (corners[4], corners[1])   # left
        ]
        
        intersections = []
        for (i, (p1, p2)) in enumerate(edges)
            f1 = f(p1[1], p1[2])
            f2 = f(p2[1], p2[2])
            
            if f1 * f2 <= 0 && !(abs(f1) < 1e-10 && abs(f2) < 1e-10)
                # Edge crosses interface
                if abs(f1) < 1e-10
                    push!(intersections, p1)
                elseif abs(f2) < 1e-10
                    push!(intersections, p2)
                else
                    # Find intersection point by linear interpolation
                    t = f1 / (f1 - f2)  # Parameter along edge
                    ix = p1[1] + t * (p2[1] - p1[1])
                    iy = p1[2] + t * (p2[2] - p1[2])
                    push!(intersections, (ix, iy))
                end
            end
        end
        
        # Construct polygon vertices
        if isempty(interior_corners)
            # No interior corners, use intersections
            if length(intersections) >= 3
                vertices = intersections
            else
                # Not enough intersections for a valid polygon
                xex[1] = (x[1] + x[2])/2
                xex[2] = (y[1] + y[2])/2
                xex[end] = 1.0
                return 0.0  # Approximation for very small cut
            end
        else
            # Combine interior corners and intersections
            vertices = vcat(interior_corners, intersections)
        end
        
        # Order vertices to form a proper polygon
        center_x = sum(v[1] for v in vertices) / length(vertices)
        center_y = sum(v[2] for v in vertices) / length(vertices)
        sort!(vertices, by = p -> atan(p[2] - center_y, p[1] - center_x))
        
        # Calculate area
        wet_area = SimpleVOF.polygon_area(vertices)
        
        # Calculate centroid using the formula
        cx, cy = SimpleVOF.polygon_centroid(vertices)

        # Calculate interface length
        interface_length = SimpleVOF.compute_interface_length(f, (x[1], x[2]), (y[1], y[2]))
        
        # Set barycenter values
        xex[1] = cx
        xex[2] = cy
        xex[end] = interface_length
        
        # Return volume
        return wet_area
    
    # 3D cases 
    # Case 6-8: 2D slices in 3D space
    # Case 6: x fixed, y and z intervals (2D slice in 3D)
    elseif length(args) == 3 && args[1] isa Number && 
           args[2] isa SVector{2} && args[3] isa SVector{2}
        x, y, z = args
        # Use SimpleVOF's 2D implementation directly
        slice_f = (yt, zt) -> f(x, yt, zt)
        frac = SimpleVOF.getcc(slice_f, (y[1], y[2]), (z[1], z[2]))
        vol = frac * (y[2] - y[1]) * (z[2] - z[1])
        
        # For cut cells, compute the accurate barycenter in the slice
        if 0 < frac < 1
            # Find the wet polygon and its centroid
            corners = [
                (y[1], z[1]),  # bottom-left
                (y[2], z[1]),  # bottom-right
                (y[2], z[2]),  # top-right
                (y[1], z[2])   # top-left
            ]
            
            # Evaluate level set at corners
            corner_values = [slice_f(yc, zc) for (yc, zc) in corners]
            
            # Find interior corners and intersections
            interior_corners = [corners[i] for i in 1:4 if corner_values[i] < 0]
            
            # Find intersections on edges
            edges = [
                (corners[1], corners[2]),  # bottom
                (corners[2], corners[3]),  # right
                (corners[3], corners[4]),  # top
                (corners[4], corners[1])   # left
            ]
            
            intersections = []
            for (p1, p2) in edges
                f1 = slice_f(p1[1], p1[2])
                f2 = slice_f(p2[1], p2[2])
                
                if f1 * f2 <= 0 && !(abs(f1) < 1e-10 && abs(f2) < 1e-10)
                    # Find intersection using linear interpolation
                    t = f1 / (f1 - f2)  # Parameter along edge
                    iy = p1[1] + t * (p2[1] - p1[1])
                    iz = p1[2] + t * (p2[2] - p1[2])
                    push!(intersections, (iy, iz))
                end
            end
            
            # Construct polygon vertices
            vertices = vcat(interior_corners, intersections)
            
            if length(vertices) >= 3
                # Calculate accurate centroid
                cy, cz = SimpleVOF.polygon_centroid(vertices)
                
                # Set barycenter values
                xex[1] = x
                xex[2] = cy
                xex[3] = cz
                xex[end] = 1.0  # Interface present
            else
                # Fallback for degenerate cases
                xex[1] = x
                xex[2] = (y[1] + y[2])/2
                xex[3] = (z[1] + z[2])/2
                xex[end] = 1.0
            end
        else
            # Full or empty cell
            xex[1] = x
            xex[2] = (y[1] + y[2])/2
            xex[3] = (z[1] + z[2])/2
            xex[end] = 0.0
        end
        
        return vol
    
    # Case 7: y fixed, x and z intervals
    elseif length(args) == 3 && args[1] isa SVector{2} && 
           args[2] isa Number && args[3] isa SVector{2}
        # Similar implementation as Case 6 but with y fixed
        x, y, z = args
        slice_f = (xt, zt) -> f(xt, y, zt)
        frac = SimpleVOF.getcc(slice_f, (x[1], x[2]), (z[1], z[2]))
        vol = frac * (x[2] - x[1]) * (z[2] - z[1])
        
        # Calculate barycenter for cut cells
        if 0 < frac < 1
            # Find the wet polygon and its centroid (similar approach)
            corners = [
                (x[1], z[1]),  # bottom-left
                (x[2], z[1]),  # bottom-right
                (x[2], z[2]),  # top-right
                (x[1], z[2])   # top-left
            ]
            
            # Evaluate level set at corners
            corner_values = [slice_f(xc, zc) for (xc, zc) in corners]
            
            # Find interior corners and intersections
            interior_corners = [corners[i] for i in 1:4 if corner_values[i] < 0]
            
            # Find intersections on edges
            edges = [
                (corners[1], corners[2]),  # bottom
                (corners[2], corners[3]),  # right
                (corners[3], corners[4]),  # top
                (corners[4], corners[1])   # left
            ]
            
            intersections = []
            for (p1, p2) in edges
                f1 = slice_f(p1[1], p1[2])
                f2 = slice_f(p2[1], p2[2])
                
                if f1 * f2 <= 0 && !(abs(f1) < 1e-10 && abs(f2) < 1e-10)
                    # Find intersection using linear interpolation
                    t = f1 / (f1 - f2)  # Parameter along edge
                    ix = p1[1] + t * (p2[1] - p1[1])
                    iz = p1[2] + t * (p2[2] - p1[2])
                    push!(intersections, (ix, iz))
                end
            end
            
            # Construct polygon vertices
            vertices = vcat(interior_corners, intersections)
            
            if length(vertices) >= 3
                # Calculate accurate centroid
                cx, cz = SimpleVOF.polygon_centroid(vertices)
                
                # Set barycenter values
                xex[1] = cx
                xex[2] = y
                xex[3] = cz
                xex[end] = 1.0  # Interface present
            else
                # Fallback for degenerate cases
                xex[1] = (x[1] + x[2])/2
                xex[2] = y
                xex[3] = (z[1] + z[2])/2
                xex[end] = 1.0
            end
        else
            # Full or empty cell
            xex[1] = (x[1] + x[2])/2
            xex[2] = y
            xex[3] = (z[1] + z[2])/2
            xex[end] = 0.0
        end
        
        return vol
    
    # Case 8: z fixed, x and y intervals
    elseif length(args) == 3 && args[1] isa SVector{2} && 
           args[2] isa SVector{2} && args[3] isa Number
        # Similar implementation as Cases 6-7 but with z fixed
        x, y, z = args
        slice_f = (xt, yt) -> f(xt, yt, z)
        frac = SimpleVOF.getcc(slice_f, (x[1], x[2]), (y[1], y[2]))
        vol = frac * (x[2] - x[1]) * (y[2] - y[1])
        
        # Calculate barycenter for cut cells
        if 0 < frac < 1
            # Find the wet polygon and its centroid (similar approach)
            corners = [
                (x[1], y[1]),  # bottom-left
                (x[2], y[1]),  # bottom-right
                (x[2], y[2]),  # top-right
                (x[1], y[2])   # top-left
            ]
            
            # Evaluate level set at corners
            corner_values = [slice_f(xc, yc) for (xc, yc) in corners]
            
            # Find interior corners and intersections
            interior_corners = [corners[i] for i in 1:4 if corner_values[i] < 0]
            
            # Find intersections on edges
            edges = [
                (corners[1], corners[2]),  # bottom
                (corners[2], corners[3]),  # right
                (corners[3], corners[4]),  # top
                (corners[4], corners[1])   # left
            ]
            
            intersections = []
            for (p1, p2) in edges
                f1 = slice_f(p1[1], p1[2])
                f2 = slice_f(p2[1], p2[2])
                
                if f1 * f2 <= 0 && !(abs(f1) < 1e-10 && abs(f2) < 1e-10)
                    # Find intersection using linear interpolation
                    t = f1 / (f1 - f2)  # Parameter along edge
                    ix = p1[1] + t * (p2[1] - p1[1])
                    iy = p1[2] + t * (p2[2] - p1[2])
                    push!(intersections, (ix, iy))
                end
            end
            
            # Construct polygon vertices
            vertices = vcat(interior_corners, intersections)
            
            if length(vertices) >= 3
                # Calculate accurate centroid
                cx, cy = SimpleVOF.polygon_centroid(vertices)
                
                # Set barycenter values
                xex[1] = cx
                xex[2] = cy
                xex[3] = z
                xex[end] = 1.0  # Interface present
            else
                # Fallback for degenerate cases
                xex[1] = (x[1] + x[2])/2
                xex[2] = (y[1] + y[2])/2
                xex[3] = z
                xex[end] = 1.0
            end
        else
            # Full or empty cell
            xex[1] = (x[1] + x[2])/2
            xex[2] = (y[1] + y[2])/2
            xex[3] = z
            xex[end] = 0.0
        end
        
        return vol
    
    # Case 9: 3D cell - Use the new accurate implementation
    elseif length(args) == 3 && all(arg isa SVector for arg in args)
        x, y, z = args
        # Use SimpleVOF's 3D implementation directly
        frac = SimpleVOF.getcc(f, (x[1], x[2]), (y[1], y[2]), (z[1], z[2]))
        vol = frac * (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])
        
        # For cut cells, compute the barycenter
        if 0 < frac < 1
            # The full barycenter calculation is already implemented in SimpleVOF
            # We need to extract it without recomputing everything
            interface_area = SimpleVOF.compute_interface_area(f, (x[1], x[2]), (y[1], y[2]), (z[1], z[2]))
            
            # Quick check of cell type
            cell_type = SimpleVOF.getcelltype(f, (x[1], x[2]), (y[1], y[2]), (z[1], z[2]))
            
            if cell_type == -1  # Cut cell
                # For now, compute an approximate barycenter
                # This is a major simplification - in practice you'd want to extract
                # the barycenter calculation from the full 3D algorithm
                
                # The 8 vertices of the cell
                corners = [
                    (x[1], y[1], z[1]),  # 1: bottom-left-back
                    (x[2], y[1], z[1]),  # 2: bottom-right-back
                    (x[2], y[2], z[1]),  # 3: top-right-back
                    (x[1], y[2], z[1]),  # 4: top-left-back
                    (x[1], y[1], z[2]),  # 5: bottom-left-front
                    (x[2], y[1], z[2]),  # 6: bottom-right-front
                    (x[2], y[2], z[2]),  # 7: top-right-front
                    (x[1], y[2], z[2])   # 8: top-left-front
                ]
                
                # Evaluate level set at corners
                corner_values = [f(cx, cy, cz) for (cx, cy, cz) in corners]
                
                # Find interior corners (ls < 0)
                interior_corners = [corners[i] for i in 1:8 if corner_values[i] < 0]
                
                if !isempty(interior_corners)
                    # If we have interior corners, use their center as an approximation
                    cx = sum(p[1] for p in interior_corners) / length(interior_corners)
                    cy = sum(p[2] for p in interior_corners) / length(interior_corners)
                    cz = sum(p[3] for p in interior_corners) / length(interior_corners)
                    
                    xex[1] = cx
                    xex[2] = cy
                    xex[3] = cz
                else
                    # No interior corners - use cell center
                    xex[1] = (x[1] + x[2])/2
                    xex[2] = (y[1] + y[2])/2
                    xex[3] = (z[1] + z[2])/2
                end
                
                xex[end] = interface_area
            else
                # Full or empty cell
                xex[1] = (x[1] + x[2])/2
                xex[2] = (y[1] + y[2])/2
                xex[3] = (z[1] + z[2])/2
                xex[end] = 0.0
            end
        else
            # Full or empty cell
            xex[1] = (x[1] + x[2])/2
            xex[2] = (y[1] + y[2])/2
            xex[3] = (z[1] + z[2])/2
            xex[end] = 0.0
        end
        
        return vol
    end
    
    # If no case matched
    error("SimpleVOF: Non implémenté pour ces arguments: $(typeof(args))")
end