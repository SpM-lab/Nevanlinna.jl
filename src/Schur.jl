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
#=
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
=#


    function calc_g(i)
        z::Complex{BigFloat} = reals.freq[i]
        param::Complex{BigFloat} = 0.0+0.0*im
        for k in 1:H
            param += ab_coeff[k]*Nevanlinna.hardy_basis(z,k-1)
            param += ab_coeff[k+H]*conj(Nevanlinna.hardy_basis(z,k-1))
        end

        theta = (abcd[1,1,i]*param + abcd[1,2,i]) / (abcd[2,1,i]*param + abcd[2,2,i])
        green = im * (one(BigFloat) + theta) / (one(BigFloat) - theta)
    end

    A = collect(imag(calc_g(i))/pi for i in 1:reals.N_real)

    second_square = 0.0
    tot_int = sum(A)*((2.0*reals.omega_max)/reals.N_real)

    #prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    #for i in 1:Int64(reals.N_real/2)
    #    prefft_spec[i] = Float64(A[Int64(reals.N_real/2)+i])
    #    prefft_spec[Int64(reals.N_real/2)+i] = Float64(A[i])
    #end

    #fft_spec = bfft(prefft_spec)*2.0*reals.omega_max/reals.N_real
    fft_spec = bfft(Float64.(A))*2.0*reals.omega_max/reals.N_real

    function calc_second(i)
        (2*pi*(i-1)/(2.0*reals.omega_max))^4 * abs(fft_spec[i])^2 /(2.0*reals.omega_max)
    end

    second_der = sum(collect(calc_second(i) for i in 2:Int64(reals.N_real/2)))

    second_der = second_der*2

    lambda::Float64 = 1e-6
    func::Complex{Float64} = abs(1-tot_int)^2 + lambda*second_der

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
