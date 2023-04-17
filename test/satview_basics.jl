import Unitful: rad, °

# These tests are the one found in the satview_basics.jl notebook
@testset "Satview Basics" begin
    # LLA
    @test LLA(10°,10°,1000) ≈ LLA((10+100*eps())*°,10°,1000)
    @test LLA(90°,10°,1000) ≈ LLA(90°,130°,1000)
    @test LLA(40°,-180°,1000) ≈ LLA(40°,180°,1000)
    @test LLA(10°,10°,1000) !== LLA((10+100*eps())*°,10°,1000)
    @test isnan(LLA(1,1,NaN))
    
    # ERA
    @test ERA(10°,1000,20°) == ERA(10°,1km,deg2rad(20)*rad)
    @test ERA(90°,1000,20°) ≈ ERA(90°,1km,deg2rad(90)*rad)
    @test isnan(ERA(1,1,NaN))


    # Geod Inverse
    em = EarthModel()
    lla1 = LLA(0°, 0°, 0)
    lla2 = LLA(1°, 0°, 0)
    @test all(geod_inverse(em.geod, lla1, lla2) .≈ (111194.92664455874, 0.0, 0.0))
end