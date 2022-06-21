@testset "util.second_deriv" begin
    # f(x) = x^2
    # f''(x) = 2
    xmax = 3
    N = 1000
    #x = collect(LinRange(0, xmax, N))
    x = collect(LinRange(0, sqrt(xmax), N)) .^ 2 # Non-uniform mesh

    res = Nevanlinna.second_deriv(x, x.^2)
    #println(x)
    #println(res)
    @test all(isapprox.(res, 2.0, rtol=0, atol=1e-5))
end

@testset "util.integrate_squared_second_deriv" begin
    # f(x) = x^3
    # f''(x) = 6*x
    # ∫_0^xmax (f''(x))^2 =12 xmax^3
    xmax = 3
    N = 10000
    #x = collect(LinRange(0, xmax, N))
    x = collect(LinRange(0, sqrt(xmax), N)) .^ 2 # Non-uniform mesh

    coeff = im

    res = Nevanlinna.integrate_squared_second_deriv(x, coeff .* x.^3)
    #println(res)
    @test isapprox(res, 12*xmax^3, atol=0, rtol=1e-3)
    #println(12 * xmax^3)
end