using TelecomUtils
using Test

@testset verbose=true "TelecomUtils.jl" begin
    include("satview_basics.jl")
    include("satview_transformations.jl")
end
