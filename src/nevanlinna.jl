function calc_phis(imags::ImagDomainData)::Vector{Complex{BigFloat}}
    phis  = Array{Complex{BigFloat}}(undef, imags.N_imag) 
    abcds = Array{Complex{BigFloat}}(undef, 2, 2, imags.N_imag) 
    phis[1] = imags.val[1]
    
    for i in 1:imags.N_imag
        view(abcds,:,:,i) .= Matrix{Complex{BigFloat}}(I, 2, 2)
    end
    
    for j in 1:imags.N_imag-1
        for k in j+1:imags.N_imag
            prod = Array{Complex{BigFloat}}(undef, 2, 2) 
            prod[1,1] = (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[2,2] = one(BigFloat)
            view(abcds,:,:,k) .= view(abcds,:,:,k)*prod
        end
        phis[j+1] = (-abcds[2,2,j+1]*imags.val[j+1]+abcds[1,2,j+1]) / (abcds[2,1,j+1]*imags.val[j+1]-abcds[1,1,j+1])
    end
    
    return phis
end

function calc_abcd(imags::ImagDomainData, 
                   reals::RealDomainData, 
                   phis::Vector{Complex{BigFloat}})::Array{Complex{BigFloat},3}
    abcd = Array{Complex{BigFloat}}(undef, 2, 2, reals.N_real) 

    for i in 1:reals.N_real
        result = Matrix{Complex{BigFloat}}(I, 2, 2) 
        z::Complex{BigFloat} = reals.freq[i]
        for j in 1:imags.N_imag
            prod = Array{Complex{BigFloat}}(undef, 2, 2)
            prod[1,1] = (z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[2,2] = one(BigFloat)
            result *= prod
        end

        abcd[:,:,i] .= result
    end
    return abcd
end

function calc_hardy_matrix(reals::RealDomainData, 
                           H::Int64)::Array{Complex{BigFloat}, 2}
    hardy_matrix = Array{Complex{BigFloat}}(undef, reals.N_real, 2*H)
    for k in 1:H
        hardy_matrix[:,k]   .=      hardy_basis.(reals.freq,k-1)
        hardy_matrix[:,k+H] .= conj(hardy_basis.(reals.freq,k-1))
    end
    return hardy_matrix
end

function calc_functional(reals::RealDomainData, 
                         abcd::Array{Complex{BigFloat},3}, 
                         H::Int64, 
                         ab_coeff::Vector{Complex{BigFloat}}, 
                         hardy_matrix::Array{Complex{BigFloat},2})::Float64
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    green = im * (one(BigFloat) .+ theta) ./ (one(BigFloat) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = sum(A)*((2.0*reals.omega_max)/reals.N_real)

    fft_spec = bfft(A)*2.0*reals.omega_max/reals.N_real

    preder_spec = fft_spec[1:Int64(reals.N_real/2)]
    t_vec = 2*pi*Vector(0:Int64(reals.N_real/2)-1)/(2.0*reals.omega_max)

    second_der = 2*sum(t_vec.^4 .* abs.(preder_spec).^2 /(2*reals.omega_max))

    lambda::Float64 = 1e-5
    func::Float64 = abs(1-tot_int)^2 + lambda*second_der

    return func
end

function evaluation(reals::RealDomainData, 
                    abcd::Array{Complex{BigFloat},3}, 
                    H::Int64, ab_coeff::Vector{Complex{BigFloat}}, 
                    hardy_matrix::Array{Complex{BigFloat},2})
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])

    reals.val .= im * (one(BigFloat) .+ theta) ./ (one(BigFloat) .- theta)
end
