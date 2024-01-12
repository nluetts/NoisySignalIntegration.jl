using NoisySignalIntegration: band_center

function band_center_testdata()
    # gaussian curve generated in LibreOffice
    X = [
        0	4.53999297624849E-05;
        1	0.000303539138078867;
        2	0.00166155727317393;
        3	0.00744658307092434;
        4	0.0273237224472926;
        5	0.0820849986238988;
        6	0.201896517994655;
        7	0.406569659740599;
        8	0.670320046035639;
        9	0.90483741803596;
        10	1;
        11	0.90483741803596;
        12	0.670320046035639;
        13	0.406569659740599;
        14	0.201896517994655;
        15	0.0820849986238988;
        16	0.0273237224472926;
        17	0.00744658307092434;
        18	0.00166155727317393;
        19	0.000303539138078867;
        20	4.53999297624849E-05;
    ]
    Curve(X[:, 1], X[:, 2])
end

@testset "gaussian peak band center" begin
crv = band_center_testdata()
    # without baseline
    @test band_center(crv, 0.0, 20.0) == 10.0
    @test band_center(crv, 6.5, 13.5) == 10.0
    @test band_center(crv, 6.5, 15.5) ≈ 10.0809916179 atol=1e-10
    @test band_center(crv, 6.3, 11.85) ≈ 9.540840963609 atol=1e-10
    # with baseline subtraction
    @test band_center(crv, 6.3, 11.85, (0.025, -0.43)) ≈ 9.3215429658485 atol=1e-10
    # this is a regression test, since I did not manually calculate the result
    @test band_center(crv, 6.3, 11.85, true) ≈ 9.354283686417842 atol=1e-10
end

@testset "4th order polynomial band center" begin
    crv = let
        f(x) = (x - 2)^4
        x = -5:0.0001:5 |> collect
        y = f.(x)
        Curve(x, y)
    end

    # result derived from wolframalpha, after baseline subtraction
    # (when integrating from -4 to 4, the endpoint to endpoint baseline
    # is 160x - 656);
    # wolframalpha query: `integrate x*(-(x-2)^4 - (160x - 656)) from -4 to 4` (numerator)
    # and:                `integrate   (-(x-2)^4 - (160x - 656)) from -4 to 4` (denominator)
    expected = (-32768/15)/(18432/5)

    @test band_center(crv, -4.0, 4.0, true) ≈ expected # true = subtract endpoint to endpoint baseline
end
