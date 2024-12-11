using ImplicitIntegration

# Utilities
isfull(val, full_val) = isapprox(val, full_val; atol=1e-8)
isempty(val) = isapprox(val, 0.0; atol=1e-8)

########################
# 1D Implementation
########################

function implicit_integration(mesh::Tuple{Vector}, Φ)
    x_coords = mesh[1]
    nx = length(x_coords)-1
    dx = x_coords[2] - x_coords[1]

    # Compute volumes (1D length segments)
    V = zeros(nx)
    for i in 1:nx
        a = (x_coords[i],)
        b = (x_coords[i+1],)
        V[i] = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
    end

    # Classify cells
    cell_types = similar(V, Int)
    for i in 1:nx
        if isempty(V[i])
            cell_types[i] = 0
        elseif isfull(V[i], dx)
            cell_types[i] = 1
        else
            cell_types[i] = -1
        end
    end

    # Compute cell centroids
    C_ω = zeros(nx) # cell centroids
    for i in 1:nx
        a = (x_coords[i],)
        b = (x_coords[i+1],)
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b).val / area
            C_ω[i] = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
        else
            C_ω[i] = 0.5*(x_coords[i]+x_coords[i+1])
        end
    end

    # Compute interface centroids
    C_γ = zeros(nx)
    Γ = zeros(nx)
    for i in 1:nx
        a = (x_coords[i],)
        b = (x_coords[i+1],)
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true).val / area
            C_γ[i] = x_c
            Γ[i] = area
        else
            C_γ[i] = NaN
            Γ[i] = 0.0
        end
    end

    # Compute W (staggered volumes)
    Wx = zeros(nx+1)
    for i in 1:nx+1
        xi = x_coords[max(i-1,1)]
        xip1 = x_coords[min(i,nx+1)]
        Wx[i] = ImplicitIntegration.integrate((x)->1, Φ, (xi,), (xip1,)).val
    end
    W = (Wx,)

    # Compute A (faces capacity)
    Ax = zeros(nx+1)
    for i in 1:nx+1
        Ax[i] = Φ((x_coords[i],)) ≤ 0.0 ? 1.0 : 0.0
    end
    A = (Ax,)

    # Compute B (values at cell centroids):
    Bx = zeros(nx)
    for i in 1:nx
        xi = C_ω[i]
        Bx[i] = Φ((xi,)) ≤ 0.0 ? 1.0 : 0.0
    end
    B = (Bx,)

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end

########################
# 2D Implementation
########################

function implicit_integration(mesh::Tuple{Vector,Vector}, Φ)
    x_coords, y_coords = mesh
    nx = length(x_coords)-1
    ny = length(y_coords)-1
    dx = x_coords[2]-x_coords[1]
    dy = y_coords[2]-y_coords[1]

    # Volumes (areas)
    V = zeros(nx, ny)
    for i in 1:nx
        for j in 1:ny
            a = (x_coords[i], y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])
            V[i,j] = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
        end
    end

    # Cell types
    cell_types = similar(V, Int)
    for i in 1:nx
        for j in 1:ny
            vol = V[i,j]
            if isempty(vol)
                cell_types[i,j] = 0
            elseif isfull(vol, dx*dy)
                cell_types[i,j] = 1
            else
                cell_types[i,j] = -1
            end
        end
    end

    # Cell centroids
    C_ω = Array{Tuple{Float64,Float64}}(undef, nx, ny)
    for i in 1:nx
        for j in 1:ny
            a = (x_coords[i], y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])
            area = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
            if area > 0
                x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b).val / area
                y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b).val / area
                x_c = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
                y_c = isnan(y_c) ? 0.5*(y_coords[j]+y_coords[j+1]) : y_c
                C_ω[i,j] = (x_c,y_c)
            else
                C_ω[i,j] = (0.5*(x_coords[i]+x_coords[i+1]), 0.5*(y_coords[j]+y_coords[j+1]))
            end
        end
    end

    # Interface centroids
    C_γ = Array{Union{Nothing,Tuple{Float64,Float64}}}(undef, nx, ny)
    Γ = zeros(nx, ny)
    for i in 1:nx
        for j in 1:ny
            a = (x_coords[i], y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])
            area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true).val
            if area > 0
                x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true).val / area
                y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b, surface=true).val / area
                C_γ[i,j] = (x_c,y_c)
                Γ[i,j] = area
            else
                C_γ[i,j] = nothing
                Γ[i,j] = 0.0
            end
        end
    end

    # Staggered volumes W
    Wx = zeros(nx+1, ny)
    Wy = zeros(nx, ny+1)
    for i in 1:nx+1
        for j in 1:ny
            xi = x_coords[max(i-1,1)]
            xip1 = x_coords[min(i,nx+1)]
            yj = y_coords[j]
            yjp1 = y_coords[j+1]
            Wx[i,j] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj), (xip1,yjp1)).val
        end
    end
    for i in 1:nx
        for j in 1:ny+1
            xi = x_coords[i]
            xip1 = x_coords[i+1]
            yj = y_coords[max(j-1,1)]
            yjp1 = y_coords[min(j,ny+1)]
            Wy[i,j] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj), (xip1,yjp1)).val
        end
    end
    W = (Wx, Wy)

    # Face capacities A (Ax, Ay)
    Ax = zeros(nx+1, ny)
    Ay = zeros(nx, ny+1)
    # For Ax (line x fixed, integrate over y):
    for i in 1:nx+1
        xi = x_coords[i]
        for j in 1:ny
            yj = y_coords[j]
            yjp1 = y_coords[j+1]
            # 1D integration in y with x fixed:
            ϕ_1d = y -> Φ((xi,y[1]))
            Ax[i,j] = ImplicitIntegration.integrate((y)->1.0, ϕ_1d, (yj,),(yjp1,)).val
        end
    end
    # For Ay (line y fixed, integrate over x):
    for j in 1:ny+1
        yj = y_coords[j]
        for i in 1:nx
            xi = x_coords[i]
            xip1 = x_coords[i+1]
            ϕ_1d = x -> Φ((x[1],yj))
            Ay[i,j] = ImplicitIntegration.integrate((x)->1.0, ϕ_1d, (xi,),(xip1,)).val
        end
    end
    A = (Ax, Ay)

    # B (Bx, By) - integrate similarly but at cell centroids or using given logic
    Bx = zeros(nx, ny)
    By = zeros(nx, ny)
    for i in 1:nx
        xi = C_ω[i,1]
        for j in 1:ny
            yj_min = y_coords[j]
            yj_max = y_coords[j+1]
            Φ_1d = (y) -> Φ((xi,y[1]))
            Bx[i,j] = ImplicitIntegration.integrate((y)->1.0, Φ_1d, (yj_min,),(yj_max,)).val
        end
    end
    for j in 1:ny
        yj = C_ω[1,j]
        for i in 1:nx
            xi_min = x_coords[i]
            xi_max = x_coords[i+1]
            Φ_1d = (x) -> Φ((x[1],yj))
            By[i,j] = ImplicitIntegration.integrate((x)->1.0, Φ_1d, (xi_min,),(xi_max,)).val
        end
    end
    B = (Bx, By)

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end

########################
# 3D Implementation
########################

function implicit_integration(mesh::Tuple{Vector,Vector,Vector}, Φ)
    x_coords, y_coords, z_coords = mesh
    nx = length(x_coords)-1
    ny = length(y_coords)-1
    nz = length(z_coords)-1
    dx = x_coords[2]-x_coords[1]
    dy = y_coords[2]-y_coords[1]
    dz = z_coords[2]-z_coords[1]

    # Volume
    V = zeros(nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        a = (x_coords[i], y_coords[j], z_coords[k])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])
        V[i,j,k] = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
    end

    # Cell types
    cell_types = similar(V, Int)
    for i in 1:nx, j in 1:ny, k in 1:nz
        vol = V[i,j,k]
        if isempty(vol)
            cell_types[i,j,k] = 0
        elseif isfull(vol, dx*dy*dz)
            cell_types[i,j,k] = 1
        else
            cell_types[i,j,k] = -1
        end
    end

    # Cell centroids
    C_ω = Array{Tuple{Float64,Float64,Float64}}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        a = (x_coords[i], y_coords[j], z_coords[k])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b).val / area
            y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b).val / area
            z_c = ImplicitIntegration.integrate((x)->x[3], Φ, a, b).val / area
            x_c = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
            y_c = isnan(y_c) ? 0.5*(y_coords[j]+y_coords[j+1]) : y_c
            z_c = isnan(z_c) ? 0.5*(z_coords[k]+z_coords[k+1]) : z_c
            C_ω[i,j,k] = (x_c,y_c,z_c)
        else
            C_ω[i,j,k] = (0.5*(x_coords[i]+x_coords[i+1]),
                          0.5*(y_coords[j]+y_coords[j+1]),
                          0.5*(z_coords[k]+z_coords[k+1]))
        end
    end

    # Interface centroids
    C_γ = Array{Union{Nothing,Tuple{Float64,Float64,Float64}},3}(undef, nx, ny, nz)
    Γ = zeros(nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        a = (x_coords[i], y_coords[j], z_coords[k])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true).val / area
            y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b, surface=true).val / area
            z_c = ImplicitIntegration.integrate((x)->x[3], Φ, a, b, surface=true).val / area
            C_γ[i,j,k] = (x_c,y_c,z_c)
            Γ[i,j,k] = area
        else
            C_γ[i,j,k] = nothing
            Γ[i,j,k] = 0.0
        end
    end

    # Staggered volumes W = (Wx, Wy, Wz)
    Wx = zeros(nx+1, ny, nz)
    Wy = zeros(nx, ny+1, nz)
    Wz = zeros(nx, ny, nz+1)
    for i in 1:nx+1, j in 1:ny, k in 1:nz
        xi = x_coords[max(i-1,1)]
        xip1 = x_coords[min(i,nx+1)]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zj = z_coords[k]
        zjp1 = z_coords[k+1]
        Wx[i,j,k] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zj), (xip1,yjp1,zjp1)).val
    end
    for i in 1:nx, j in 1:ny+1, k in 1:nz
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[max(j-1,1)]
        yjp1 = y_coords[min(j,ny+1)]
        zj = z_coords[k]
        zjp1 = z_coords[k+1]
        Wy[i,j,k] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zj), (xip1,yjp1,zjp1)).val
    end
    for i in 1:nx, j in 1:ny, k in 1:nz+1
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zj = z_coords[max(k-1,1)]
        zjp1 = z_coords[min(k,nz+1)]
        Wz[i,j,k] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zj), (xip1,yjp1,zjp1)).val
    end
    W = (Wx, Wy, Wz)

    # A = (Ax, Ay, Az), integrating over faces with one coordinate fixed
    # Here we just replicate logic similar to 2D, extended to 3D
    Ax = zeros(nx+1, ny, nz)
    for i in 1:nx+1
        for j in 1:ny 
            for k in 1:nz
                xi = x_coords[i]
                y0, y1 = y_coords[j], y_coords[j+1]
                z0, z1 = z_coords[k], z_coords[k+1]
                ϕ_2d = (y,z) -> Φ((xi,y,z))
                Ax[i,j,k] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (y0,z0), (y1,z1)).val
            end
        end
    end

    Ay = zeros(nx, ny+1, nz)
    for i in 1:nx
        for j in 1:ny+1
            for k in 1:nz
                yj = y_coords[j]
                x0, x1 = x_coords[i], x_coords[i+1]
                z0, z1 = z_coords[k], z_coords[k+1]
                ϕ_2d = (x,z) -> Φ((x,yj,z))
                Ay[i,j,k] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (x0,z0), (x1,z1)).val
            end
        end
    end

    Az = zeros(nx, ny, nz+1)
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz+1
                zk = z_coords[k]
                x0, x1 = x_coords[i], x_coords[i+1]
                y0, y1 = y_coords[j], y_coords[j+1]
                ϕ_2d = (x,y) -> Φ((x,y,zk))
                Az[i,j,k] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (x0,y0), (x1,y1)).val
            end
        end
    end
    A = (Ax, Ay, Az)

    # B = (Bx, By, Bz)
    # For simplicity, evaluate Φ at cell centroid and check sign
    Bx = zeros(nx, ny, nz)
    By = zeros(nx, ny, nz)
    Bz = zeros(nx, ny, nz)
    for i in 1:nx
        xi = C_ω[i,1,1]
        for j in 1:ny
            yj = y_coords[j]
            yjmax = y_coords[j+1]
            for k in 1:nz
                zk = z_coords[k]
                zkmax = z_coords[k+1]
                Φ_2d = (y,z) -> Φ((xi,y,z))
                Bx[i,j,k] = ImplicitIntegration.integrate((y,z)->1.0, Φ_2d, (yj,zk), (yjmax,zkmax)).val
            end
        end
    end

    for j in 1:ny
        yj = C_ω[1,j,1]
        for i in 1:nx
            xi = x_coords[i]
            ximax = x_coords[i+1]
            for k in 1:nz
                zk = z_coords[k]
                zkmax = z_coords[k+1]
                Φ_2d = (x,z) -> Φ((x,yj,z))
                By[i,j,k] = ImplicitIntegration.integrate((x,z)->1.0, Φ_2d, (xi,zk), (ximax,zkmax)).val
            end
        end
    end

    for k in 1:nz
        zk = C_ω[1,1,k]
        for i in 1:nx
            xi = x_coords[i]
            ximax = x_coords[i+1]
            for j in 1:ny
                yj = y_coords[j]
                yjmax = y_coords[j+1]
                Φ_2d = (x,y) -> Φ((x,y,zk))
                Bz[i,j,k] = ImplicitIntegration.integrate((x,y)->1.0, Φ_2d, (xi,yj), (ximax,yjmax)).val
            end
        end
    end

    B = (Bx, By, Bz)

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end
