using CartesianGeometry
using Plots
Plots.default(show = true)
const T = Float64

# Définition de la grille
xyz = (0:0.1:2, 0:0.1:2)

# Définition de la fonction levelset
levelset = (x, y, _=0) -> x^2 + y^2 - 1

# Plot
plot(xyz..., levelset, st = [:contourf], fill = true, c = :viridis)
readline()

# Étape 1 : Calcul initial des volumes et des barycentres
V, bary = integrate(Tuple{0}, levelset, xyz, T, nan)

# Restructurer les données
dims = ntuple(i -> length(xyz[i]) , length(xyz))
V = reshape(V, dims)
bary = reshape(bary, dims)
@show size(V)
# Étape 2 : Définition du seuil
threshold = 0.001

# Étape 3 : Trouver les petites cellules
small_cells = find_small_cells(V, threshold)
@show size(small_cells)
# Étape 4 : Fusion des petites cellules avec les voisines
for ind in small_cells
    neighbors = find_neighbors(ind, dims)
    # Filtrer les voisins valides (volume > 0)
    valid_neighbors = filter(n -> V[n] > 0, neighbors)
    if !isempty(valid_neighbors)
        # Sélectionner le voisin avec le plus grand volume
        neighbor_ind = valid_neighbors[argmax(V[valid_neighbors])]
        merge_cells!(V, bary, ind, neighbor_ind)
    end
end

# Étape 5 : Recalculer As, Ws, Bs avec les nouvelles données

# Recalculer As (adapter la fonction integrate pour qu'elle retourne des matrices)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

# Restructurer As si nécessaire
# ...

# Recalculer Ws et Bs en utilisant les nouvelles données
Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

@show size(V)