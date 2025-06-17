using CartesianGeometry
using Test
using Statistics
using Printf
using LinearAlgebra

"""
Compare VOFI and SimpleVOF methods for volume and surface fractions in 3D
"""
function compare_methods_3d()
    println("Comparing VOFI vs SimpleVOF methods in 3D")
    println("=========================================")
    
    # Generate 3D grid with increasing resolution
    resolutions = [11, 21]  # Lower resolution for 3D to keep runtime reasonable
    
    # Test case: Sphere centered at (0.5, 0.5, 0.5) with radius 0.3
    center = (0.5, 0.5, 0.5)
    radius = 0.3
    levelset = HyperSphere(radius, center)
    
    # Analytical values for comparison
    analytical_volume = (4/3) * π * radius^3
    analytical_surface = 4π * radius^2
    
    for resolution in resolutions
        println("\n## Resolution: $resolution x $resolution x $resolution")
        
        x = collect(range(0, 1, length=resolution))
        y = collect(range(0, 1, length=resolution))
        z = collect(range(0, 1, length=resolution))
        xyz = (x, y, z)
        
        # Run VOFI method
        println("\nRunning VOFI method...")
        
        # First kind - volume fractions
        V_vofi, bary_vofi, interface_vofi, types_vofi = 
            integrate(Tuple{0}, levelset, xyz, Float64, nan; method=:vofi)
        
        # First kind - surface fractions
        As_vofi = integrate(Tuple{1}, levelset, xyz, Float64, nan; method=:vofi)
        
        # Second kind - volume fractions using barycenters
        Ws_vofi = integrate(Tuple{0}, levelset, xyz, Float64, nan, bary_vofi; method=:vofi)
        
        # Second kind - surface fractions using barycenters
        Bs_vofi = integrate(Tuple{1}, levelset, xyz, Float64, nan, bary_vofi; method=:vofi)
        
        # Run SimpleVOF method
        println("\nRunning SimpleVOF method...")
        
        # First kind - volume fractions
        V_simple, bary_simple, interface_simple, types_simple = 
            integrate(Tuple{0}, levelset, xyz, Float64, nan; method=:simple)
        
        # First kind - surface fractions
        As_simple = integrate(Tuple{1}, levelset, xyz, Float64, nan; method=:simple)
        
        # Second kind - volume fractions using barycenters
        Ws_simple = integrate(Tuple{0}, levelset, xyz, Float64, nan, bary_simple; method=:simple)
        
        # Second kind - surface fractions using barycenters
        Bs_simple = integrate(Tuple{1}, levelset, xyz, Float64, nan, bary_simple; method=:simple)
        
        # Calculate total volume and surface area
        vofi_volume = sum(filter(!isnan, V_vofi))
        simple_volume = sum(filter(!isnan, V_simple))

        println(V_simple)
        println(V_vofi)
        println(V_vofi - V_simple)
        
        # Function to print difference vectors with limited output
        function print_diff_vector(diff_vector, name, max_display=10)
            println("\n  $name difference vector:")
            nonzero_diff = diff_vector[abs.(diff_vector) .> 1e-6]
            if length(nonzero_diff) > 0
                if length(nonzero_diff) <= max_display
                    println("    ", round.(nonzero_diff, digits=6))
                else
                    println("    First $max_display elements: ", round.(nonzero_diff[1:max_display], digits=6))
                    println("    ... ($(length(nonzero_diff)-max_display) more elements)")
                end
                println("    Number of non-zero differences: $(length(nonzero_diff))")
            else
                println("    No significant differences found")
            end
        end
        
        # Compare volume fractions (V)
        valid_indices = .!isnan.(V_vofi)
        volume_diff = V_vofi[valid_indices] - V_simple[valid_indices]
        type_diff = types_vofi[valid_indices] - types_simple[valid_indices]
        
        println("\n1. VOLUME FRACTIONS COMPARISON:")
        println("-------------------------------")
        println("  Analytical volume:    ", @sprintf("%.6f", analytical_volume))
        println("  VOFI total volume:    ", @sprintf("%.6f", vofi_volume))
        println("  SimpleVOF volume:     ", @sprintf("%.6f", simple_volume))
        println("  VOFI error:           ", @sprintf("%.6f", abs(vofi_volume - analytical_volume)))
        println("  SimpleVOF error:      ", @sprintf("%.6f", abs(simple_volume - analytical_volume)))
        println("  Cell-wise differences:")
        println("    Mean difference:    ", @sprintf("%.6f", mean(abs.(volume_diff))))
        println("    Max difference:     ", @sprintf("%.6f", maximum(abs.(volume_diff))))
        println("    Different types:    ", count(type_diff .!= 0), " of ", length(type_diff))
        
        # Print V difference vector
        print_diff_vector(volume_diff, "V")
        
        # Compare surface fractions (As)
        println("\n2. SURFACE FRACTIONS COMPARISON:")
        println("---------------------------------")
        for i in 1:length(As_vofi)
            valid_indices = .!isnan.(As_vofi[i])
            area_diff = As_vofi[i][valid_indices] - As_simple[i][valid_indices]
            
            vofi_surface = sum(filter(!isnan, As_vofi[i]))
            simple_surface = sum(filter(!isnan, As_simple[i]))
            
            println("  Direction $i:")
            println("    VOFI total surface:   ", @sprintf("%.6f", vofi_surface))
            println("    SimpleVOF surface:    ", @sprintf("%.6f", simple_surface))
            println("    Mean difference:      ", @sprintf("%.6f", mean(abs.(area_diff))))
            println("    Max difference:       ", @sprintf("%.6f", maximum(abs.(area_diff))))
            
            # Print A difference vector
            print_diff_vector(area_diff, "A[$i]")
        end
        
        # Calculate total surface area from all directions
        vofi_total_surface = sum(sum(filter(!isnan, As_vofi[i])) for i in 1:length(As_vofi))
        simple_total_surface = sum(sum(filter(!isnan, As_simple[i])) for i in 1:length(As_simple))
        
        println("\n  Total Surface Area Summary:")
        println("    Analytical surface:   ", @sprintf("%.6f", analytical_surface))
        println("    VOFI total surface:   ", @sprintf("%.6f", vofi_total_surface))
        println("    SimpleVOF surface:    ", @sprintf("%.6f", simple_total_surface))
        println("    VOFI error:           ", @sprintf("%.6f", abs(vofi_total_surface - analytical_surface)))
        println("    SimpleVOF error:      ", @sprintf("%.6f", abs(simple_total_surface - analytical_surface)))
        
        # Compare second kind volume fractions (Ws)
        println("\n3. SECOND KIND VOLUME FRACTIONS:")
        println("--------------------------------")
        for i in 1:length(Ws_vofi)
            valid_indices = .!isnan.(Ws_vofi[i])
            ws_diff = Ws_vofi[i][valid_indices] - Ws_simple[i][valid_indices]
            
            println("  Direction $i:")
            println("    Mean difference:      ", @sprintf("%.6f", mean(abs.(ws_diff))))
            println("    Max difference:       ", @sprintf("%.6f", maximum(abs.(ws_diff))))
            
            # Print W difference vector
            print_diff_vector(ws_diff, "W[$i]")
        end
        
        # Compare second kind surface fractions (Bs)
        println("\n4. SECOND KIND SURFACE FRACTIONS:")
        println("---------------------------------")
        for i in 1:length(Bs_vofi)
            valid_indices = .!isnan.(Bs_vofi[i])
            bs_diff = Bs_vofi[i][valid_indices] - Bs_simple[i][valid_indices]
            
            println("  Direction $i:")
            println("    Mean difference:      ", @sprintf("%.6f", mean(abs.(bs_diff))))
            println("    Max difference:       ", @sprintf("%.6f", maximum(abs.(bs_diff))))
            
            # Print B difference vector
            print_diff_vector(bs_diff, "B[$i]")
        end
        """
        # Compare barycenters
        println("\n5. BARYCENTERS COMPARISON:")
        println("--------------------------")
        for i in axes(bary_vofi, 1)
            # For SVectors, we need to check if any component is NaN
            valid_indices = findall(x -> !any(isnan, x), bary_vofi[i,:])
            
            if !isempty(valid_indices)
                # Get the values as vectors for comparison
                vofi_values = [bary_vofi[i,j] for j in valid_indices]
                simple_values = [bary_simple[i,j] for j in valid_indices]
                
                # Calculate component-wise differences for each dimension
                differences = [norm(vofi - simple) for (vofi, simple) in zip(vofi_values, simple_values)]
                
                println("  Dimension $i:")
                println("    Mean difference:      ", @sprintf("%.6f", mean(differences)))
                println("    Max difference:       ", @sprintf("%.6f", maximum(differences)))
                println("    Num comparisons:      ", length(differences))
            else
                println("  Dimension $i: No valid comparisons")
            end
        end
        """
    end
end

# Run the test
compare_methods_3d()