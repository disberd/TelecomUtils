import Unitful: rad, °

# These tests are the one found in the satview_basics.jl notebook
@testset "Satview Basics" begin
    # LLA
    @test LLA(10°,10°,1000) ≈ LLA((10+100*eps())*°,10°,1000)
    @test LLA(90°,10°,1000) ≈ LLA(90°,130°,1000)
    @test LLA(40°,-180°,1000) ≈ LLA(40°,180°,1000)
    @test LLA(0°, 0°, 0km) ≉ LLA(1.1e-5°, 0°, 0km)
    @test LLA(0°, 0°, 0km) ≈ LLA(1e-5°, 0°, 0km)
    @test LLA(0°, 0°, 0km) ≈ LLA(1e-5°, 1e-5°, 1e-6km)
    @test LLA(0°, 0°, 0km) ≉ LLA(1e-5°, 1e-5°, 1.1e-6km)
    @test_throws "atol" isapprox(LLA(0°, 0°, 0km),LLA(1e-5°, 1e-5°, 1e-6km); atol = 0.2)
    @test LLA(10°,10°,1000) !== LLA((10+100*eps())*°,10°,1000)
    @test isnan(LLA(1,1,NaN))
    
    # ERA
    @test ERA(10°,1000,20°) == ERA(10°,1km,deg2rad(20)*rad)
    @test ERA(90°,1000,20°) ≈ ERA(90°,1km,deg2rad(90)*rad)
    @test isnan(ERA(1,1,NaN))
    @test ERA(0°, 0km, 0°) ≉ ERA(1.1e-5°, 0km, 0°)
    @test ERA(0°, 0km, 0°) ≈ ERA(1e-5°, 0km, 0°)
    @test ERA(0°, 0km, 0°) ≈ ERA(1e-5°, 0km + 1e-6km, 0°)
    @test ERA(0°, 0km, 0°) ≉ ERA(1e-5°, 0km + 1.1e-6km, 0°)
    @test_throws "atol" isapprox(ERA(0°, 0km, 0°),ERA(1e-5°, 0km, 0°); atol = 0.2)


    # Geod Inverse
    em = EarthModel()
    lla1 = LLA(0°, 0°, 0)
    lla2 = LLA(1°, 0°, 0)
    @test all(geod_inverse(em.geod, lla1, lla2) .≈ (111194.92664455874, 0.0, 0.0))
end