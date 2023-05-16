import Unitful: °, km
using Rotations
import LinearAlgebra: normalize
import TelecomUtils: earth_intersection, ExtraOutput, to_radians
using TelecomUtils
using StaticArrays
using Test

sp_ell = SphericalEllipsoid()
### ECEF <-> LLA ###
lla2ecef = ECEFfromLLA(sp_ell)
@test lla2ecef(LLA(0, 0, 0)) ≈ [6371e3, 0, 0]
lla1 = LLA(30°, 45°, 100km)
@test lla2ecef(lla1) |> inv(lla2ecef) ≈ lla1 # Test transformation inversion

### Earth Intersection ###
@testset "Earth Intersection" begin
    sat_lla = LLA(0, 0, 500km)
    sat_ecef = lla2ecef(sat_lla)

    n̂ = normalize(-sat_ecef)
    intersection_ecef = earth_intersection(n̂, sat_ecef, sp_ell.a, sp_ell.b)
    @test intersection_ecef ≈ [sp_ell.a, 0, 0]
    # Test that if the direction is opposite of the earth the result is NaN
    @test all(map(isnan, earth_intersection(-n̂, sat_ecef, sp_ell.a, sp_ell.b)))
end

### Compute sat positions ###
@testset "Compute Sat Positions" begin
    h = 6371e3 + 735e3
    lla = LLA(10°, 25°, 0km)
    era2ecef = ECEFfromERA(lla)
    el, az = 35°, 80°
    satecef = compute_sat_position(era2ecef, el, az; h)
    invera = inv(era2ecef)(satecef)
    @test invera.el ≈ el && invera.az ≈ az
end

@testset "UVfromLLA" begin
    sat_lla = LLA(0°, 0°, 600km)
    target_lla = [
        LLA(0°, 1°, 0km), # Right - U Negative, V 0
        LLA(1°, 0°, 0km), # Top - U 0, V Positive
        LLA(0°, -1°, 0km), # Left - U Positive, V 0
        LLA(-1°, 0°, 0km), # Bottom - U 0, V Negative
    ]
    target_uv = map(target_lla) do lla
        UVfromLLA(sat_lla; ellipsoid=SphericalEllipsoid())(lla) |> normalize
    end
    @test target_uv[1] == [-1, 0]
    @test target_uv[2] == [0, 1]
    @test target_uv[3] == [1, 0]
    @test target_uv[4] == [0, -1]


    # We now test that targets behind the reference direction are not visible (NaN)
    sat_lla = LLA(0, 0, 600km)
    target_lla = LLA(0, 0, 610km)
    lla2uv = UVfromLLA(sat_lla)
    target_uv = lla2uv(target_lla)
    @test all(isnan.(target_uv))
    # We test that if we rotate by 180 degrees around Y the matrix we actually find a valid uv
    lla2uv = UVfromLLA(lla2uv.origin, RotY(180°) * lla2uv.R, lla2uv.ellipsoid)
    @test all(map(!isnan, lla2uv(target_lla)))
end

@testset "LLAfromUV" begin
    sp = SphericalEllipsoid()
    l2e = ECEFfromLLA(sp)
    # We find the pointing that corresponds to Edge of Earth and check various combination in its vicinity
    sat_lla = LLA(0, 0, 600km)
    sat_ecef = l2e(sat_lla)
    eoe_scan = asin(sp.a / (sp.a + sat_lla.alt))
    u = sin(eoe_scan)
    uv2lla = LLAfromUV(sat_lla; ellipsoid=sp)
    @test !isnan(uv2lla((u * (1 - eps()), 0))) # We should find a solution because we are pointing slightly less than EoE
    @test isnan(uv2lla((u * (1 + eps()), 0))) # We should not find a solution because we are pointing slightly more than EoE
    @test !isnan(uv2lla((u * (1 - eps()), 0), h=100e3)) # We should find a solution because we are looking at 100km above earth
    @test !isnan(uv2lla((u * (1 + eps()), 0), h=100e3)) # We should find a solution because we are looking at 100km above earth
    @test isnan(uv2lla((u * (1 - eps()), 0), h=700e3)) # We should not find a solution because we are looking at 100km above the satellite alitude and with an angle slightly lower than eoe scan, so the corresponding valid point in the pointing direction is located behind earth
    @test !isnan(uv2lla((u * (1 + eps()), 0), h=700e3)) # We should find a solution because we are pointing more than eoe_scan so the earth is not blocking the view of the corresponding point


    lla2uv = inv(uv2lla)
    target_uv = SA_F64[0.1, 0.1]
    target_lla, r = uv2lla(target_uv, ExtraOutput())
    uv2, r2 = lla2uv(target_lla, ExtraOutput())
    @test uv2 ≈ target_uv
    @test r2 ≈ r
end

# Add angle offset
@testset "Add Angle Offset" begin
    # Test adding having an offset on the same direction
    θ_vec = range(10°, 80°; step=20°)
    φ_vec = range(0°, 359°; step=30°)
    function test_φ(x, y)
        x̂ = rem2pi(to_radians(x), RoundNearest)
        ŷ = rem2pi(to_radians(y), RoundNearest)
        result = abs(x̂) ≈ abs(ŷ) ≈ π || Base.isapprox(x̂, ŷ; atol=1e-10, rtol=1e-5)
        result || @info "Phi" x y x̂ ŷ Base.isapprox(x̂, ŷ; atol=1e-10, rtol=1e-5)
        result
    end
    for θ in θ_vec
        for φ in φ_vec
            p_add = add_angle_offset((θ, φ), (5°, φ); input_type=:thetaphi, output_type=:thetaphi)
            p_sub = add_angle_offset((θ, φ), (5°, φ + 180°); input_type=:thetaphi, output_type=:thetaphi)
            @test p_add[1] ≈ θ + 5° && test_φ(p_add[2], φ)
            @test p_sub[1] ≈ θ - 5° && test_φ(p_sub[2], φ)
        end
    end

    # Test a perpendicular offset
    for θ₁ in range(10°, 50°; step=20°)
        for θ₂ in range(10°, 50°; step=20°)
            for φ in range(0°, 270°; step=45°)
                p = add_angle_offset((θ₁, φ - 45°), (θ₂, φ + 45°); input_type=:thetaphi, output_type=:thetaphi)
                # Since they are perpendicular, we can  use the right spherical triangle rule (cosine law)
                @test p[1] ≈ acos(cos(θ₁) * cos(θ₂))
            end
        end
    end

    # Test that pointing behind throws an error
    @test_throws "behind the viewer" add_angle_offset((0.7, 0), (50°, 0))
end
