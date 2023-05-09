using TelecomUtils
using Test
using StaticArrays

@testset verbose=true "TelecomUtils.jl" begin
    Base.isnan(v::StaticArray) = any(isnan, v)
    include("refview_basics.jl")
    include("refview_transformations.jl")
    include("refview_struct.jl")
end
