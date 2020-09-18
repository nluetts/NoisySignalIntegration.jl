@testset "allapproxequal()" begin
    @test nsi.allapproxequal([])                     == true
    @test nsi.allapproxequal([1.0])                  == true
    @test nsi.allapproxequal([1.00000001, 1.0, 1.0]) == true
    # tolerance is 1e-8          1234567
    @test nsi.allapproxequal([ 1.0000001, 1.0, 1.0]) == false
    @test nsi.allapproxequal([ 1, 1, 1])             == true
    @test nsi.allapproxequal([ 2, 1, 1])             == false
    @test_throws MethodError nsi.allapproxequal([ "A", "A"])
end


@testset "verify_same_length()" begin
    @test nsi.verify_same_length([], [])                 === nothing
    @test nsi.verify_same_length(["A"], ["B"])           === nothing
    @test nsi.verify_same_length(["A", "B"], ["C", "D"]) === nothing
    @test nsi.verify_same_length([1], [2])               === nothing
    @test nsi.verify_same_length([1, 2], [3, 4])         === nothing
    @test nsi.verify_same_length([1 2; 3 4], [5 6; 7 8]) === nothing
    @test_throws ArgumentError nsi.verify_same_length([1], [])                
    @test_throws ArgumentError nsi.verify_same_length(["A"], [])          
    @test_throws ArgumentError nsi.verify_same_length(["A", "B"], ["C"])
    @test_throws ArgumentError nsi.verify_same_length([1], [])              
    @test_throws ArgumentError nsi.verify_same_length([1, 2], [3])
    @test_throws ArgumentError nsi.verify_same_length([1 2; 3 4], [5 6;])
end


@testset "left_right_from_peak()" begin
    x = collect(Float64, 1:20)
    y = zeros(Float64, 20)
    y[10] = 1
    @test nsi.left_right_from_peak(x, y, 10, 1) == [9.5, 10.5]
    @test nsi.left_right_from_peak(x, y, 10, 0.7) == [9.65, 10.35]
    y[10] = 0
    y[1]  = 1
    @test nsi.left_right_from_peak(x, y, 1, 1) == [0.5, 1.5]
    y[1] = 0
    y[20]  = 1
    @test nsi.left_right_from_peak(x, y, 20, 1) == [19.5, 20.5]
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
            @test nsi.trapz(x, f.(x), l, r, subtract_baseline=false) ≈ ref_1 rtol=1e-7
            @test nsi.trapz(x, f.(x), l, r) ≈ ref_2 rtol=1e-7
        end
    end

    begin
        x = collect(Float64, 1:10)
        y = collect(Float64, 1:10)
        @test nsi.trapz(x, y, 1, 3, subtract_baseline=false) == 4.0
        @test nsi.trapz(x, y, 1.1, 2.1, subtract_baseline=false) ≈ 1.6
        @test nsi.trapz(x, y, 1.1, 1.2, subtract_baseline=false) ≈ 1.1*0.1 + 0.1*0.1*0.5
        @test nsi.trapz(x, y, 1.0, 1.1, subtract_baseline=false) ≈ 0.1 + 0.1*0.1*0.5
        @test nsi.trapz(x, y, 3.9, 4.0, subtract_baseline=false) ≈ 3.9*0.1 + 0.1*0.1*0.5
        @test nsi.trapz(x, y, 0, 10) == 0.0
        @test nsi.trapz(x, y, 1, 11) == 0.0
        @test nsi.trapz(x, y, 0, 11) == 0.0
    end

    begin
        x = [1, 2, 3, 4, 6, 7.0]
        y = [1, 2, 4, 2.5, 3, 2]
        @test nsi.trapz(x, y, 1, 7, subtract_baseline=false) == 15.75
        @test nsi.trapz(x, y, 0, 10, subtract_baseline=false) == 15.75
        @test nsi.trapz(x, y, 1, 2, subtract_baseline=false) == 1.5
        @test nsi.trapz(x, y, 3, 4, subtract_baseline=false) == 3.25
        @test nsi.trapz(x, y, 1, 6, subtract_baseline=false) == 13.25
        @test nsi.trapz(x, y, 6, 1, subtract_baseline=false) == 13.25
        @test nsi.trapz(x, y, 1.5, 7, subtract_baseline=false) == 15.125
    end
end

@testset "@samples" begin
    seed!(1)
    @testset "Normal" begin
        ref = Particles(10_000, Normal(1, 1))
        out = @samples 10_000 1 ± 1
        @test mean(ref) ≈ mean(out)
        @test std(ref) ≈ std(out)
        # pass expression to macro
        out = @samples 10_000 (2 / 2) ± (3 - 2)
        @test mean(ref) ≈ mean(out)
        @test std(ref) ≈ std(out)
        # pass variables to macro
        out = let
            a = 1
            b = 1
            @samples 10_000 a ± b
        end
        @test mean(ref) ≈ mean(out)
        @test std(ref) ≈ std(out)
    end
    @testset "Uniform" begin
        ref = Particles(10_000, Uniform(1, 2))
        out = @samples 10_000 1 .. 2
        @test mean(ref) ≈ mean(out)
        @test std(ref) ≈ std(out)
        # pass expression to macro
        out = @samples 10_000 (2 / 2) .. (3 - 1)
        @test mean(ref) ≈ mean(out)
        @test std(ref) ≈ std(out)
        # pass variables to macro
        out = let
            a = 1
            b = 2
            @samples 10_000 a .. b
        end
        @test mean(ref) ≈ mean(out)
        @test std(ref) ≈ std(out)
    end
    @test_throws ArgumentError (@samples 1 1 + 1)
    @test_throws ArgumentError (@samples 1 1 + 1 + 1)
end