@testset "moment" begin
 
    T = BigFloat
    setprecision(2048)

    #define spectral function
    gaussian(x, mu, sigma) = exp(-0.5*((x-mu)/sigma)^2)/(sqrt(2*π)*sigma)
    rho(omega) = 0.5*gaussian(omega, 2.0, 1.0) + 0.5*gaussian(omega, -2.0, 1.0)


    function generate_input_data(rho::Function, beta::Float64)
        lambda = 1e+4
        wmax = lambda/beta
        basis = SparseIR.FiniteTempBasisSet(beta, wmax, 1e-15)

        rhol = [overlap(basis.basis_f.v[l], rho) for l in 1:length(basis.basis_f)]
        gl = - basis.basis_f.s .* rhol
        gw = evaluate(basis.smpl_wn_f, gl)

        hnw = length(basis.smpl_wn_f.sampling_points)÷2

        #To exclude effect of enviroment, we limit data until 30th

        input_smpl = Array{Complex{T}}(undef, 31)
        input_gw   = Array{Complex{T}}(undef, 31)

        for i in 1:31
            input_smpl[i]= SparseIR.valueim(basis.smpl_wn_f.sampling_points[hnw+i], beta)
            input_gw[i]  = gw[hnw+i]
        end

        return input_smpl, input_gw
    end

    beta = 100. #inverse temperature
    input_smpl, input_gw = generate_input_data(rho, beta)

    N_real    = 1000  #demension of array of output
    omega_max = 10.0  #energy cutoff of real axis
    eta       = 0.001 #broaden parameter 
    sum_rule  = 1.0   #sum rule
    H_max     = 50    #cutoff of Hardy basis
    lambda    = 1e-4  #regularization parameter
    iter_tol  = 1000  #upper bound of iteration
    moments = Complex{T}.([1, 0, 5, 0, 43])

    sol = Nevanlinna.HamburgerNevanlinnaSolver(moments, input_smpl, input_gw, N_real, omega_max, eta, sum_rule, H_max, iter_tol, lambda, optimization=false)

    calc_moment = Vector{Float64}(undef,5)
    for i in 1:5
        xk = real.(sol.nev_st.reals.freq).^(i-1)
        y = Float64.(imag.(sol.val))/pi .* xk
        calc_moment[i] = Float64(Nevanlinna.integrate(real.(sol.nev_st.reals.freq), y))
    end

    test_moment = ([0.9999328252706802, 2.1117052430510528e-10, 5.005359475447759, -8.294576805137473e-9, 43.293142029878446])
    @test isapprox(calc_moment, test_moment; atol = 1e-4)

    println(calc_moment)

end
