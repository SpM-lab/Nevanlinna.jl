using LinearAlgebra
using FFTW

struct ImagDomainData
    N_imag::Int64
    freq  ::Array{Complex{BigFloat},1}
    val   ::Array{Complex{BigFloat},1}
end

function ImagDomainData(N_imag::Int64,
                        matsu ::Array{BigFloat,1},
                        green ::Array{Complex{BigFloat},1}
                        )::ImagDomainData
    val  = Array{Complex{BigFloat}}(undef, N_imag) 
    freq = Array{Complex{BigFloat}}(undef, N_imag) 
    
    for i in 1:N_imag
        freq[i] = matsu[i]*im
        val[i]  = (-green[i] - im) / (-green[i] + im)
    end
    
    Pick = Array{Complex{BigFloat}}(undef, N_imag, N_imag)
    
    for j in 1:N_imag
        for i in 1:N_imag
            freq_i = (freq[i] - im) / (freq[i] + im)
            freq_j = (freq[j] - im) / (freq[j] + im)
            nom = one(BigFloat) - val[i] * conj(val[j])
            den = one(BigFloat) - freq_i * conj(freq_j)
            Pick[i,j] = nom / den
        end
        Pick[j,j] += big(1e-250)
    end
    
    success = issuccess(cholesky(Pick,check = false))
    
    if success
        println("Pick matrix is positive semi-definite.")
    else
        println("Pick matrix is non positive semi-definite matrix in Schur method.")
    end
    
    freq = reverse(freq)
    val  = reverse(val)
    
    return ImagDomainData(N_imag, freq, val)
end

function calc_phis(imag::ImagDomainData)::Vector{Complex{BigFloat}}
    phis  = Array{Complex{BigFloat}}(undef, imag.N_imag) 
    abcds = Array{Complex{BigFloat}}(undef, 2, 2, imag.N_imag) 
    phis[1] = imag.val[1]
    
    for i in 1:imag.N_imag
        view(abcds,:,:,i) .= Matrix{Complex{BigFloat}}(I, 2, 2)
    end
    
    for j in 1:imag.N_imag-1
        for k in j+1:imag.N_imag
            prod = Array{Complex{BigFloat}}(undef, 2, 2) 
            prod[1,1] = (imag.freq[k] - imag.freq[j]) / (imag.freq[k] - conj(imag.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(imag.freq[k] - imag.freq[j]) / (imag.freq[k] - conj(imag.freq[j]))
            prod[2,2] = one(BigFloat)
            view(abcds,:,:,k) .= view(abcds,:,:,k)*prod
        end
        phis[j+1] = (-abcds[2,2,j+1]*imag.val[j+1]+abcds[1,2,j+1]) / (abcds[2,1,j+1]*imag.val[j+1]-abcds[1,1,j+1])
    end
    
    return phis
end

struct RealDomainData
    N_real   ::Int64
    omega_min::Float64
    omega_max::Float64
    eta      ::Float64
    freq     ::Array{Complex{BigFloat},1}
    val      ::Array{Complex{BigFloat},1}
end

function RealDomainData(N_real   ::Int64,
                        omega_min::Float64,
                        omega_max::Float64,
                        eta      ::Float64
                        )::RealDomainData
    val  = Array{Complex{BigFloat}}(undef, N_real) 
    freq = Array{Complex{BigFloat}}(undef, N_real) 
    
    inter::BigFloat = parse(BigFloat, string((omega_max - omega_min) / (N_real-1)))
    temp ::BigFloat = parse(BigFloat, string(omega_min))
    
    freq[1] = parse(BigFloat, string(omega_min)) + parse(BigFloat, string(eta))*im
    for i in 2:N_real
        temp += inter
        freq[i] = temp + parse(BigFloat, string(eta))*im
    end
    
    return RealDomainData(N_real, omega_min, omega_max, eta, freq, val)
end

function calc_abcd(imag::ImagDomainData, reals::RealDomainData, phis::Vector{Complex{BigFloat}})::Array{Complex{BigFloat},3}
    abcd = Array{Complex{BigFloat}}(undef, 2, 2, reals.N_real) 

    for i in 1:reals.N_real
        result = Matrix{Complex{BigFloat}}(I, 2, 2) 
        z::Complex{BigFloat} = reals.freq[i]
        for j in 1:imag.N_imag
            prod = Array{Complex{BigFloat}}(undef, 2, 2)
            prod[1,1] = (z - imag.freq[j]) / (z - conj(imag.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imag.freq[j]) / (z - conj(imag.freq[j]))
            prod[2,2] = one(BigFloat)
            result *= prod
        end

        abcd[:,:,i] .= result
    end
    return abcd
end

function calc_functional(reals::RealDomainData, abcd::Array{Complex{BigFloat},3}, H::Int64, ab_coeff::Vector{Float64} )
    function calc_g(i)
        z::Complex{BigFloat} = reals.freq[i]
        param::Complex{BigFloat} = 0.0+0.0*im
        for k in 1:H
            param += ab_coeff[k]*hardy_basis(z,k-1)
            param += ab_coeff[k+H]*conj(hardy_basis(z,k-1))
        end
        
        theta = (abcd[1,1,i]*param + abcd[1,2,i]) / (abcd[2,1,i]*param + abcd[2,2,i])
        green = im * (one(BigFloat) + theta) / (one(BigFloat) - theta)
    end

    A = collect(imag(calc_g(i))/pi for i in 1:reals.N_real)

    lambda::Float64 = 1e-6
    tot_int::Float64 = 0.0
    second_square::Float64 = 0.0
    for i in 1:reals.N_real
        tot_int += A[i]*((reals.omega_max-reals.omega_min)/reals.N_real)
        second_der_i = (A[mod1(i+1,reals.N_real)] + A[mod1(i-1,reals.N_real)] - 2*A[i])*((reals.omega_max-reals.omega_min)/reals.N_real)^(-2)
        second_square += (second_der_i)^2*((reals.omega_max-reals.omega_min)/reals.N_real)
    end
    
    func::Float64 = (1-tot_int) + lambda*second_square
    

end



function evaluation(imag::ImagDomainData, reals::RealDomainData, phis::Vector{Complex{BigFloat}}, H::Int64, ab_coeff::Vector{Float64} )
    for i in 1:reals.N_real
        result = Matrix{Complex{BigFloat}}(I, 2, 2) 
        z::Complex{BigFloat} = reals.freq[i]
        for j in 1:imag.N_imag
            prod = Array{Complex{BigFloat}}(undef, 2, 2)
            prod[1,1] = (z - imag.freq[j]) / (z - conj(imag.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imag.freq[j]) / (z - conj(imag.freq[j]))
            prod[2,2] = one(BigFloat)
            result *= prod
        end
        
        param::Complex{BigFloat} = 0.0+0.0*im
        for k in 1:H
            param += ab_coeff[k]*hardy_basis(z,k-1)
            param += ab_coeff[k+H]*conj(hardy_basis(z,k-1))
        end
        
        theta = (result[1,1]*param + result[1,2]) / (result[2,1]*param + result[2,2])
        reals.val[i] = im * (one(BigFloat) + theta) / (one(BigFloat) - theta)
    end
end

function hardy_basis(z::Complex{BigFloat}, k::Int64)
    (z-im)^k/(sqrt(pi)*(z+im)^(k+1))
end

function hardy_basis(x::Float64, y::Float64, k::Int64)
    bigx = parse(BigFloat, string(x))
    bigy = parse(BigFloat, string(y))
    
    z::Complex{BigFloat} = x +im*y
    
    return hardy_basis(z, k)
end

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
