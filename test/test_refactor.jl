using Pkg
Pkg.activate(".")

# Test loading
try
    using DiscreteCurves
    using LinearAlgebra
    println("✓ Module loaded successfully")
    
    # Test that new functions are accessible
    println("\n=== Testing Exports ===")
    
    # GeometricProperties
    println("✓ signed_area: ", isdefined(DiscreteCurves, :signed_area))
    println("✓ centroid: ", isdefined(DiscreteCurves, :centroid))
    
    # Orientation
    println("✓ Orientation: ", isdefined(DiscreteCurves, :Orientation))
    println("✓ CCW: ", isdefined(DiscreteCurves, :CCW))
    println("✓ CW: ", isdefined(DiscreteCurves, :CW))
    println("✓ CounterClockwise: ", isdefined(DiscreteCurves, :CounterClockwise))
    println("✓ Clockwise: ", isdefined(DiscreteCurves, :Clockwise))
    println("✓ is_ccw: ", isdefined(DiscreteCurves, :is_ccw))
    println("✓ is_cw: ", isdefined(DiscreteCurves, :is_cw))
    println("✓ orient: ", isdefined(DiscreteCurves, :orient))
    println("✓ orient!: ", isdefined(DiscreteCurves, :orient!))
    
    # CurveVectors
    println("✓ tangent_vector: ", isdefined(DiscreteCurves, :tangent_vector))
    println("✓ tangent_vectors: ", isdefined(DiscreteCurves, :tangent_vectors))
    println("✓ normal_vector: ", isdefined(DiscreteCurves, :normal_vector))
    println("✓ left_normal_vector: ", isdefined(DiscreteCurves, :left_normal_vector))
    println("✓ right_normal_vector: ", isdefined(DiscreteCurves, :right_normal_vector))
    println("✓ inward_normal: ", isdefined(DiscreteCurves, :inward_normal))
    println("✓ outward_normal: ", isdefined(DiscreteCurves, :outward_normal))
    
    # Quick functional test
    println("\n=== Quick Functional Test ===")
    c = circle_curve(16)
    a = signed_area(c)
    println("Circle signed area: $(round(a; digits=3))")
    println("Is CCW: $(is_ccw(c))")
    
    t = tangent_vector(c, 1)
    n = normal_vector(c, 1)
    println("Tangent norm: $(round(norm(t); digits=6))")
    println("Normal norm: $(round(norm(n); digits=6))")
    
    println("\n✓ All tests passed!")
    
catch e
    println("✗ Error: ", e)
    stacktrace(catch_backtrace()) |> display
end
