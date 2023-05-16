using SafeTestsets

@safetestset "RefView Basics" begin include("refview_basics.jl") end
@safetestset "RefView Transformations" begin include("refview_transformations.jl") end
@safetestset "RefView Struct" begin include("refview_struct.jl") end
