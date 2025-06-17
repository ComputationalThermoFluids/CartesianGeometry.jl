using CartesianGeometry
using Test
using Statistics
using Printf

"""
Compare VOFI and SimpleVOF methods for volume fractions in 1D
"""
function compare_methods_1d()
    println("Comparing VOFI vs SimpleVOF methods in 1D")
    println("=========================================")
    
    # Generate grids with increasing resolution
    resolutions = [21, 51, 101]
    
    # Test case 1: Interval [-1,1] with interface at x=0
    println("\n## TEST CASE 1: Linear function with zero at x=0")
    
    # Simple level set: f(x) = x
    levelset_linear = x -> x
    
    # Analytical values for comparison
    analytical_volume_1 = 1.0  # Length of interval on negative side
    
    for resolution in resolutions
        println("\n### Resolution: $resolution points")
        
        # Create 1D grid
        x = collect(range(-1, 1, length=resolution))
        
        # Run VOFI method
        println("\nRunning VOFI method...")
        V_vofi, bary_vofi, interface_vofi, types_vofi = 
            integrate(Tuple{0}, levelset_linear, (x,), Float64, nan; method=:vofi)
        
        # Run SimpleVOF method
        println("Running SimpleVOF method...")
        V_simple, bary_simple, interface_simple, types_simple = 
            integrate(Tuple{0}, levelset_linear, (x,), Float64, nan; method=:simple)
        
        # Calculate total volume
        vofi_volume = sum(filter(!isnan, V_vofi))
        simple_volume = sum(filter(!isnan, V_simple))
        
        # Compare volume fractions
        valid_indices = .!isnan.(V_vofi)
        volume_diff = V_vofi[valid_indices] - V_simple[valid_indices]
        type_diff = types_vofi[valid_indices] - types_simple[valid_indices]
        
        println("\n1. VOLUME FRACTIONS COMPARISON:")
        println("  Analytical volume:    ", @sprintf("%.6f", analytical_volume_1))
        println("  VOFI total volume:    ", @sprintf("%.6f", vofi_volume))
        println("  SimpleVOF volume:     ", @sprintf("%.6f", simple_volume))
        println("  VOFI error:           ", @sprintf("%.6f", abs(vofi_volume - analytical_volume_1)))
        println("  SimpleVOF error:      ", @sprintf("%.6f", abs(simple_volume - analytical_volume_1)))
        println("  Cell-wise differences:")
        println("    Mean difference:    ", @sprintf("%.6f", mean(abs.(volume_diff))))
        println("    Max difference:     ", @sprintf("%.6f", maximum(abs.(volume_diff))))
        println("    Different types:    ", count(type_diff .!= 0), " of ", length(type_diff))
        
        # Compare interface positions
        println("\n2. INTERFACE POSITIONS COMPARISON:")
  
    end
    
    # Test case 2: Cubic function with multiple interfaces
    println("\n## TEST CASE 2: Cubic function with multiple interfaces")
    
    # Level set: f(x) = x(x-0.5)(x+0.5)
    levelset_cubic = x -> x * (x - 0.5) * (x + 0.5)
    
    # Analytical values (roots are at x = -0.5, 0, 0.5)
    # Negative regions are [-0.5, 0] and [0.5, 1]
    analytical_volume_2 = 0.5 + 0.5  # Total length in intervals: [-1,-0.5] and [0,0.5]
    
    for resolution in resolutions
        println("\n### Resolution: $resolution points")
        
        # Create 1D grid
        x = collect(range(-1, 1, length=resolution))
        
        # Run VOFI method
        println("\nRunning VOFI method...")
        V_vofi, bary_vofi, interface_vofi, types_vofi = 
            integrate(Tuple{0}, levelset_cubic, (x,), Float64, nan; method=:vofi)
        
        # Run SimpleVOF method
        println("Running SimpleVOF method...")
        V_simple, bary_simple, interface_simple, types_simple = 
            integrate(Tuple{0}, levelset_cubic, (x,), Float64, nan; method=:simple)

        # Calculate total volume
        vofi_volume = sum(filter(!isnan, V_vofi))
        simple_volume = sum(filter(!isnan, V_simple))
        
        # Compare volume fractions
        valid_indices = .!isnan.(V_vofi)
        volume_diff = V_vofi[valid_indices] - V_simple[valid_indices]
        type_diff = types_vofi[valid_indices] - types_simple[valid_indices]
        
        println("\n1. VOLUME FRACTIONS COMPARISON:")
        println("  Analytical volume:    ", @sprintf("%.6f", analytical_volume_2))
        println("  VOFI total volume:    ", @sprintf("%.6f", vofi_volume))
        println("  SimpleVOF volume:     ", @sprintf("%.6f", simple_volume))
        println("  VOFI error:           ", @sprintf("%.6f", abs(vofi_volume - analytical_volume_2)))
        println("  SimpleVOF error:      ", @sprintf("%.6f", abs(simple_volume - analytical_volume_2)))
        println("  Cell-wise differences:")
        println("    Mean difference:    ", @sprintf("%.6f", mean(abs.(volume_diff))))
        println("    Max difference:     ", @sprintf("%.6f", maximum(abs.(volume_diff))))
        println("    Different types:    ", count(type_diff .!= 0), " of ", length(type_diff))
        
        # Compare interface positions
        println("\n2. INTERFACE COUNTS:")
    end
    
    # Test case 3: Oscillating function (sine wave)
    println("\n## TEST CASE 3: Oscillating function (sine wave)")
    
    # Level set: f(x) = sin(2π*x)
    levelset_sine = x -> sin(2π*x)
    
    # Analytical values (negative when sin(2π*x) < 0)
    # Negative in intervals [0.25, 0.75]
    analytical_volume_3 = 0.5  # Length of negative region
    
    for resolution in resolutions
        println("\n### Resolution: $resolution points")
        
        # Create 1D grid
        x = collect(range(0, 1, length=resolution))
        
        # Run VOFI method
        println("\nRunning VOFI method...")
        V_vofi, bary_vofi, interface_vofi, types_vofi = 
            integrate(Tuple{0}, levelset_sine, (x,), Float64, nan; method=:vofi)
        
        # Run SimpleVOF method
        println("Running SimpleVOF method...")
        V_simple, bary_simple, interface_simple, types_simple = 
            integrate(Tuple{0}, levelset_sine, (x,), Float64, nan; method=:simple)
        
        # Calculate total volume
        vofi_volume = sum(filter(!isnan, V_vofi))
        simple_volume = sum(filter(!isnan, V_simple))
        
        # Compare volume fractions
        valid_indices = .!isnan.(V_vofi)
        volume_diff = V_vofi[valid_indices] - V_simple[valid_indices]
        type_diff = types_vofi[valid_indices] - types_simple[valid_indices]
        
        println("\n1. VOLUME FRACTIONS COMPARISON:")
        println("  Analytical volume:    ", @sprintf("%.6f", analytical_volume_3))
        println("  VOFI total volume:    ", @sprintf("%.6f", vofi_volume))
        println("  SimpleVOF volume:     ", @sprintf("%.6f", simple_volume))
        println("  VOFI error:           ", @sprintf("%.6f", abs(vofi_volume - analytical_volume_3)))
        println("  SimpleVOF error:      ", @sprintf("%.6f", abs(simple_volume - analytical_volume_3)))
        println("  Cell-wise differences:")
        println("    Mean difference:    ", @sprintf("%.6f", mean(abs.(volume_diff))))
        println("    Max difference:     ", @sprintf("%.6f", maximum(abs.(volume_diff))))
        println("    Different types:    ", count(type_diff .!= 0), " of ", length(type_diff))
    end
    
    # Test case 4: Sharp sigmoid function
    println("\n## TEST CASE 4: Sharp sigmoid function")
    
    # Level set: sigmoid centered at x=0.3
    sigmoid_center = 0.3
    sigmoid_steepness = 50.0  # Higher = sharper transition
    levelset_sigmoid = x -> tanh(sigmoid_steepness * (x - sigmoid_center))
    
    # Analytical values (negative when x < sigmoid_center)
    analytical_volume_4 = sigmoid_center  # Length of negative region [0, sigmoid_center]
    
    for resolution in resolutions
        println("\n### Resolution: $resolution points")
        
        # Create 1D grid
        x = collect(range(0, 1, length=resolution))
        
        # Run VOFI method
        println("\nRunning VOFI method...")
        V_vofi, bary_vofi, interface_vofi, types_vofi = 
            integrate(Tuple{0}, levelset_sigmoid, (x,), Float64, nan; method=:vofi)
        
        # Run SimpleVOF method
        println("Running SimpleVOF method...")
        V_simple, bary_simple, interface_simple, types_simple = 
            integrate(Tuple{0}, levelset_sigmoid, (x,), Float64, nan; method=:simple)
        
        # Calculate total volume
        vofi_volume = sum(filter(!isnan, V_vofi))
        simple_volume = sum(filter(!isnan, V_simple))
        
        # Compare volume fractions
        valid_indices = .!isnan.(V_vofi)
        volume_diff = V_vofi[valid_indices] - V_simple[valid_indices]
        type_diff = types_vofi[valid_indices] - types_simple[valid_indices]
        
        println("\n1. VOLUME FRACTIONS COMPARISON:")
        println("  Analytical volume:    ", @sprintf("%.6f", analytical_volume_4))
        println("  VOFI total volume:    ", @sprintf("%.6f", vofi_volume))
        println("  SimpleVOF volume:     ", @sprintf("%.6f", simple_volume))
        println("  VOFI error:           ", @sprintf("%.6f", abs(vofi_volume - analytical_volume_4)))
        println("  SimpleVOF error:      ", @sprintf("%.6f", abs(simple_volume - analytical_volume_4)))
        println("  Cell-wise differences:")
        println("    Mean difference:    ", @sprintf("%.6f", mean(abs.(volume_diff))))
        println("    Max difference:     ", @sprintf("%.6f", maximum(abs.(volume_diff))))
        println("    Different types:    ", count(type_diff .!= 0), " of ", length(type_diff))
        
       
        
    end
    
    println("\n=========================================")
    println("1D comparison tests complete!")
end

# Run the test
compare_methods_1d()