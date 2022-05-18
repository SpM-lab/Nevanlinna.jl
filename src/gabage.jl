function calc_functional(reals::RealDomainData)
    #=
    prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    for i in 1:Int64(reals.N_real/2)
        prefft_spec[i] = Float64(imag(reals.val[Int64(reals.N_real/2)+i])/pi)
        prefft_spec[Int64(reals.N_real/2)+i] = Float64(imag(reals.val[i])/pi)
    end

    fft_spec = bfft(prefft_spec)*((reals.omega_max-reals.omega_min)/reals.N_real)
    second_der::Complex{Float64} = 0.0 + 0.0*im

    for i in 2:Int64(reals.N_real/2)
        t = 2*pi*(i-1)/(reals.omega_max-reals.omega_min)
        second_der += t^4 * fft_spec[i] * fft_spec[reals.N_real-i+2] /(reals.omega_max-reals.omega_min)
    end
    
    second_der = second_der*2

    lambda::Float64 = 1e-6
    func::Complex{Float64} = abs(1-fft_spec[1])^2 + lambda*second_der
    =#
    
    lambda::Float64 = 1e-6
    tot_int::Float64 = 0.0
    second_square::Float64 = 0.0
    for i in 1:reals.N_real
        tot_int += imag(reals.val[i]/pi)*((reals.omega_max-reals.omega_min)/reals.N_real)
        second_der_i = (imag(reals.val[mod1(i+1,reals.N_real)]/pi) + imag(reals.val[mod1(i-1,reals.N_real)]/pi) - 2*imag(reals.val[i]/pi))*((reals.omega_max-reals.omega_min)/reals.N_real)^(-2)
        second_square += (second_der_i)^2*((reals.omega_max-reals.omega_min)/reals.N_real)
    end
    
    func::Float64 = (1-tot_int) + lambda*second_square
    
    #=
    prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    for i in 1:Int64(reals.N_real/2)
        prefft_spec[i] = Float64(imag(reals.val[Int64(reals.N_real/2)+i])/pi)
        prefft_spec[Int64(reals.N_real/2)+i] = Float64(imag(reals.val[i])/pi)
    end

    fft_spec = bfft(prefft_spec)*2.0*reals.omega_max/reals.N_real

    second_der_prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    for i in 1:Int64(reals.N_real/2)
        t = pi*(i-1)/reals.omega_max
        second_der_prefft_spec[i] = fft_spec[i]*t^2
        second_der_prefft_spec[reals.N_real-i+1] = fft_spec[reals.N_real-i+1]*t^2
    end

    second_der_fft_spec = fft(second_der_prefft_spec)/(2.0*reals.omega_max)

    second_der::Complex{Float64} = 0.0 + 0.0*im

    for i in 1:reals.N_real
        #second_der += second_der_fft_spec[i]*second_der_fft_spec[i]*2.0*reals.omega_max/reals.N_real
        second_der += abs(second_der_fft_spec[i])*2.0*reals.omega_max/reals.N_real
    end
    
    lambda::Float64 = 1e-4
    func::Float64 = abs(1-fft_spec[1])^2 + lambda*second_der
    =#
end
