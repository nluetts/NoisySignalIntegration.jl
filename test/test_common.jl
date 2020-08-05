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


@testset "trapz()" begin
    # Test integration function trapz using analytic integrals;
    # integral functions (F) of functions f from
    # https://en.wikipedia.org/wiki/Lists_of_integrals

    bounds = [(-1, 4), (1, 4), (-1.4, 7.9), (11.6, 17.4)]

    test_data = [
        (f = x->  exp(3x), F = x->        1/3*exp(3x)),
        (f = x->  1.2^(x), F = x->   1.2^(x)/log(1.2)),
        (f = x->   sin(x), F = x->            -cos(x)),
        (f = x-> 1/(2x+3), F = x-> 1/2*log(abs(2x+3)))
    ]
    x = collect(-20:0.0001:20)
    for t in test_data
        for (l, r) in bounds
            f, F = t
            ref_1 = F(r) - F(l)
            ref_2 = ref_1 - 1/2*(f(r) + f(l))*(r - l) # 1/2*(f(r) + f(l))*(r - l) term: subtracting baseline
            @test mci.trapz(x, f.(x), l, r, subtract_baseline=false) ≈ ref_1 rtol=1e-7
            @test mci.trapz(x, f.(x), l, r) ≈ ref_2 rtol=1e-7
        end
    end

    begin
        x = collect(Float64, 1:10)
        y = collect(Float64, 1:10)
        @test mci.trapz(x, y, 1, 3, subtract_baseline=false) == 4.0
        @test mci.trapz(x, y, 1.1, 2.1, subtract_baseline=false) ≈ 1.6
        @test mci.trapz(x, y, 1.1, 1.2, subtract_baseline=false) ≈ 1.1*0.1 + 0.1*0.1*0.5
        @test mci.trapz(x, y, 1.0, 1.1, subtract_baseline=false) ≈ 0.1 + 0.1*0.1*0.5
        @test mci.trapz(x, y, 3.9, 4.0, subtract_baseline=false) ≈ 3.9*0.1 + 0.1*0.1*0.5
        @test mci.trapz(x, y, 0, 10) == 0.0
        @test mci.trapz(x, y, 1, 11) == 0.0
        @test mci.trapz(x, y, 0, 11) == 0.0
    end

    begin
        x = [1, 2, 3, 4, 6, 7.0]
        y = [1, 2, 4, 2.5, 3, 2]
        @test mci.trapz(x, y, 1, 7, subtract_baseline=false) == 15.75
        @test mci.trapz(x, y, 0, 10, subtract_baseline=false) == 15.75
        @test mci.trapz(x, y, 1, 2, subtract_baseline=false) == 1.5
        @test mci.trapz(x, y, 3, 4, subtract_baseline=false) == 3.25
        @test mci.trapz(x, y, 1, 6, subtract_baseline=false) == 13.25
        @test mci.trapz(x, y, 6, 1, subtract_baseline=false) == 13.25
        @test mci.trapz(x, y, 1.5, 7, subtract_baseline=false) == 15.125
    end


end
