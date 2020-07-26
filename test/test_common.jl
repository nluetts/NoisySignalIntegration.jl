using Test
using MCIntegrate

const mci = MCIntegrate


@testset "allapproxequal()" begin
    @test mci.allapproxequal([])                     == true
    @test mci.allapproxequal([1.0])                  == true
    @test mci.allapproxequal([1.00000001, 1.0, 1.0]) == true
    # tolerance is 1e-8          1234567
    @test mci.allapproxequal([ 1.0000001, 1.0, 1.0]) == false
    @test mci.allapproxequal([ 1, 1, 1])             == true
    @test mci.allapproxequal([ 2, 1, 1])             == false
    @test_throws MethodError mci.allapproxequal([ "A", "A"])
end


@testset "verify_same_length()" begin
    @test mci.verify_same_length([], [])                 === nothing
    @test mci.verify_same_length(["A"], ["B"])           === nothing
    @test mci.verify_same_length(["A", "B"], ["C", "D"]) === nothing
    @test mci.verify_same_length([1], [2])               === nothing
    @test mci.verify_same_length([1, 2], [3, 4])         === nothing
    @test mci.verify_same_length([1 2; 3 4], [5 6; 7 8]) === nothing
    @test_throws ArgumentError mci.verify_same_length([1], [])                
    @test_throws ArgumentError mci.verify_same_length(["A"], [])          
    @test_throws ArgumentError mci.verify_same_length(["A", "B"], ["C"])
    @test_throws ArgumentError mci.verify_same_length([1], [])              
    @test_throws ArgumentError mci.verify_same_length([1, 2], [3])
    @test_throws ArgumentError mci.verify_same_length([1 2; 3 4], [5 6;])
end


@testset "detrend()" begin
    for n ∈ 2:4
        @test begin
            x = collect(Float64, 1:20)
            y = x.^n
            all(yᵢ + 1.0 ≈ 1.0 for yᵢ ∈ mci.detrend(x, y, n)) == true
        end
    end
    @test begin
        x = collect(Float64, 1:20)
        y = @. 0.1 + 2x - 0.01x^2
        all(yᵢ + 1.0 ≈ 1.0 for yᵢ ∈ mci.detrend(x, y, 2)) == true
    end
end


@testset "left_right_from_peak()" begin
    x = collect(Float64, 1:20)
    y = zeros(Float64, 20)
    y[10] = 1
    @test mci.left_right_from_peak(x, y, 10, 1) == [9.5, 10.5]
    @test mci.left_right_from_peak(x, y, 10, 0.7) == [9.65, 10.35]
    y[10] = 0
    y[1]  = 1
    @test mci.left_right_from_peak(x, y, 1, 1) == [0.5, 1.5]
    y[1] = 0
    y[20]  = 1
    @test mci.left_right_from_peak(x, y, 20, 1) == [19.5, 20.5]
end


@testset "get_linear_param()" begin
    @test mci.get_linear_param(0, 0.5, 1, 1.5) == (1.0, 1.0)
    @test mci.get_linear_param(3, 2, 3, 5) == (0.0, 1.0)
end


@testset "trapz()" begin
    # Test integration function trapz using analytic integrals;
    # integral functions (F) of functions f from
    # https://en.wikipedia.org/wiki/Lists_of_integrals
    test_data = [
        (f = x->  exp(3x), F = x->        1/3*exp(3x), l =    -1, r =    3),
        (f = x->  1.2^(x), F = x->   1.2^(x)/log(1.2), l =     1, r =    3),
        (f = x->   sin(x), F = x->            -cos(x), l =  -2.5, r =  7.9),
        (f = x-> 1/(2x+3), F = x-> 1/2*log(abs(2x+3)), l = -11.6, r = -2.4)
    ]
    x = collect(-20:0.0001:20)
    for t in test_data
        f, F, l, r = t
        trapz_result = mci.trapz(x, f.(x), l, r)
        # 1/2*(f(r) + f(l))*(r - l) term: subtracting baseline
        reference = F(r) - F(l) - 1/2*(f(r) + f(l))*(r - l)
        @test trapz_result ≈ reference rtol=1e-7
    end

    begin
        x = collect(Float64, 1:10)
        y = ones(Float64, length(x))
        @test mci.trapz(x, y, 1, 3) == 0.0
        @test mci.trapz(x, y, 1.1, 2.1) == 0.0
        @test mci.trapz(x, y, 1.1, 1.2) == 0.0
        @test mci.trapz(x, y, 1, 10) == 0.0
        @test mci.trapz(x, y, 0, 10) == 0.0
        @test mci.trapz(x, y, 1, 11) == 0.0
        @test mci.trapz(x, y, 0, 11) == 0.0
    end


end
