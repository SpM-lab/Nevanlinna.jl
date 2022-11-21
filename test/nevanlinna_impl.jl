@testset "nevanlinna_impl" begin
  
    T = BigFloat
    gaussian(x, mu, sigma) = exp(-((x-mu)/sigma)^2)/(sqrt(π)*sigma)

    rho(omega) = 0.8*gaussian(omega, -1.0, 1.6) + 0.2*gaussian(omega, 3, 1)
    setprecision(256)

    hnw = 38

    test_gw   = Array{Complex{T}}(undef, hnw)
    test_smpl = Array{Complex{T}}(undef, hnw)

    f = open((@__DIR__) * "/c++/result/green.dat", "r")
    for i in 1:hnw
        list = readline(f)
        s  = split(list,'\t')
        o  = parse(BigFloat, s[1])
        re = parse(BigFloat, s[2])
        ii = parse(BigFloat, s[3])
        test_smpl[i] = im*o
        test_gw[i]   = re + ii*im
    end
    close(f)

    N_real    =  6000
    omega_max =  10.0
    eta       =  0.001
    H_max     =  50
    ab_coeff  = zeros(ComplexF64, 2*H_max)
    lambda    = 1e-5
    iter_tol  = 1000
    N_imag    =  Nevanlinna.calc_opt_N_imag(hnw, test_smpl, test_gw)


    imaginary = Nevanlinna.ImagDomainData(test_smpl, test_gw, N_imag)
    raw_reals = Nevanlinna.RealDomainData(N_real, omega_max, eta, 1.0, T=T, mesh=:test)

    phis = Nevanlinna.calc_phis(imaginary)
    abcd = Nevanlinna.calc_abcd(imaginary, raw_reals, phis)
    hardy_matrix = Nevanlinna.calc_hardy_matrix(raw_reals, H_max)

    Nevanlinna.evaluation!(raw_reals, abcd, H_max, ab_coeff, hardy_matrix)

    spec = imag.(raw_reals.val)/pi

    cpp_phis = Array{Complex{T}}(undef, N_imag)
    f = open((@__DIR__) * "/c++/result/phis.dat", "r")
    for i in 1:N_imag
        list = readline(f)
        s  = split(list,'\t')
        real_phi = parse(BigFloat, s[1])
        imag_phi = parse(BigFloat, s[2])
        cpp_phis[i]   = real_phi + imag_phi*im
    end
    close(f)

    @test cpp_phis ≈ phis

    cpp_abcd = Array{Complex{T}}(undef, 2, 2, N_real)
    f = open((@__DIR__) * "/c++/result/abcd.dat", "r")
    for i in 1:N_real
        list = readline(f)
        s  = split(list,'\t')
        real_11 = parse(BigFloat, s[1])
        imag_11 = parse(BigFloat, s[2])
        real_12 = parse(BigFloat, s[3])
        imag_12 = parse(BigFloat, s[4])
        real_21 = parse(BigFloat, s[5])
        imag_21 = parse(BigFloat, s[6])
        real_22 = parse(BigFloat, s[7])
        imag_22 = parse(BigFloat, s[8])
        cpp_abcd[1,1,i]   = real_11 + imag_11*im
        cpp_abcd[1,2,i]   = real_12 + imag_12*im
        cpp_abcd[2,1,i]   = real_21 + imag_21*im
        cpp_abcd[2,2,i]   = real_22 + imag_22*im
    end
    close(f)
    
    @test cpp_abcd ≈ abcd

    cpp_spec = Array{Float64}(undef, N_real)
    f = open((@__DIR__) * "/c++/result/out_spec.dat", "r")
    for i in 1:N_real
        list = readline(f)
        s  = split(list,'\t')
        spectral = parse(Float64, s[2])
        cpp_spec[i] = spectral
    end
    close(f)

    @test (cpp_spec) ≈ Float64.(spec)
end
