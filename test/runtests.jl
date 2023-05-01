using TelecomUtils
using Test
using StaticArrays

@testset verbose=true "TelecomUtils.jl" begin
    Base.isnan(v::StaticArray) = any(isnan, v)
    include("satview_basics.jl")
    include("satview_transformations.jl")
    include("satview_struct.jl")
end
