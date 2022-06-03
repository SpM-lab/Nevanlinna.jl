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

function calc_hardy_matrix(reals::RealDomainData, H::Int64)::Array{Complex{BigFloat}, 2}
    hardy_matrix = Array{Complex{BigFloat}}(undef, reals.N_real, 2*H)
    for k in 1:H
        hardy_matrix[:,k]   .=      hardy_basis.(reals.freq,k-1)
        hardy_matrix[:,k+H] .= conj(hardy_basis.(reals.freq,k-1))
    end
    return hardy_matrix
end

function calc_functional(reals::RealDomainData, abcd::Array{Complex{BigFloat},3}, H::Int64, ab_coeff::Vector{Complex{BigFloat}}, hardy_matrix::Array{Complex{BigFloat},2})::Float64
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    green = im * (one(BigFloat) .+ theta) ./ (one(BigFloat) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = sum(A)*((2.0*reals.omega_max)/reals.N_real)

    #本来はこうすべきだが、以下のようにしても結果は変わらない
    #prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    #for i in 1:Int64(reals.N_real/2)
    #    prefft_spec[i] = Float64(A[Int64(reals.N_real/2)+i])
    #    prefft_spec[Int64(reals.N_real/2)+i] = Float64(A[i])
    #end
    #fft_spec = bfft(prefft_spec)*2.0*reals.omega_max/reals.N_real
 
    fft_spec = bfft(A)*2.0*reals.omega_max/reals.N_real

    preder_spec = fft_spec[1:Int64(reals.N_real/2)]
    t_vec = 2*pi*Vector(0:Int64(reals.N_real/2)-1)/(2.0*reals.omega_max)

    second_der = 2*sum(t_vec.^4 .* abs.(preder_spec).^2 /(2*reals.omega_max))

    lambda::Float64 = 1e-5
    func::Float64 = abs(1-tot_int)^2 + lambda*second_der

    #println(abs(1-tot_int)^2 )
    #println(second_der)

    return func
end

function evaluation(imag::ImagDomainData, reals::RealDomainData, phis::Vector{Complex{BigFloat}}, H::Int64, ab_coeff::Vector{Complex{BigFloat}})
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

function evaluation(imag::ImagDomainData, reals::RealDomainData, abcd::Array{Complex{BigFloat},3}, H::Int64, ab_coeff::Vector{Complex{BigFloat}}, hardy_matrix::Array{Complex{BigFloat},2})
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])

    reals.val .= im * (one(BigFloat) .+ theta) ./ (one(BigFloat) .- theta)
end
