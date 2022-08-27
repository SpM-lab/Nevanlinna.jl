function calc_phis(imags::ImagDomainData{T})::Vector{Complex{T}} where {T<:Real}
    phis  = Array{Complex{T}}(undef, imags.N_imag) 
    abcds = Array{Complex{T}}(undef, 2, 2, imags.N_imag) 
    phis[1] = imags.val[1]
    
    for i in 1:imags.N_imag
        view(abcds,:,:,i) .= Matrix{Complex{T}}(I, 2, 2)
    end
    
    for j in 1:imags.N_imag-1
        for k in j+1:imags.N_imag
            prod = Array{Complex{T}}(undef, 2, 2) 
            prod[1,1] = (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j]) * (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[2,2] = one(T)
            view(abcds,:,:,k) .= view(abcds,:,:,k)*prod
        end
        phis[j+1] = (-abcds[2,2,j+1]*imags.val[j+1] + abcds[1,2,j+1]) / (abcds[2,1,j+1]*imags.val[j+1] - abcds[1,1,j+1])
    end
    
    return phis
end

function calc_abcd(imags::ImagDomainData{T}, 
                   reals::RealDomainData{T}, 
                   phis::Vector{Complex{T}}
                   )::Array{Complex{T},3} where {T<:Real}
    abcd = Array{Complex{T}}(undef, 2, 2, reals.N_real) 

    for i in 1:reals.N_real
        result = Matrix{Complex{T}}(I, 2, 2) 
        z::Complex{T} = reals.freq[i]
        for j in 1:imags.N_imag
            prod = Array{Complex{T}}(undef, 2, 2)
            prod[1,1] = (z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[2,2] = one(T)
            result *= prod
        end

        abcd[:,:,i] .= result
    end
    return abcd
end

function calc_hardy_matrix(reals::RealDomainData{T}, 
                           H::Int64
                           )::Array{Complex{T}, 2} where {T<:Real}
    hardy_matrix = Array{Complex{T}}(undef, reals.N_real, 2*H)
    for k in 1:H
        hardy_matrix[:,2*k-1] .=      hardy_basis.(reals.freq,k-1)
        hardy_matrix[:,2*k]   .= conj(hardy_basis.(reals.freq,k-1))
    end
    return hardy_matrix
end

function calc_functional(reals::RealDomainData{T}, 
                         abcd::Array{Complex{T},3}, 
                         H::Int64, 
                         ab_coeff::Vector{Complex{S}}, 
                         hardy_matrix::Array{Complex{T},2};
                         lambda::Float64 = 1e-5
                         )::Float64 where {S<:Real, T<:Real}
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    green = im * (one(T) .+ theta) ./ (one(T) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = integrate(reals.freq, A)
    second_der = integrate_squared_second_deriv(reals.freq, A) 

    max_theta = findmax(abs.(param))[1]
    func = abs(reals.sum-tot_int)^2 + lambda*second_der

    return func
end

function evaluation!(reals::RealDomainData{T}, 
                    abcd::Array{Complex{T},3}, 
                    H::Int64, ab_coeff::Vector{Complex{S}}, 
                    hardy_matrix::Array{Complex{T},2}
                    )::Bool where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        println("max_theta=",max_theta)
        println("hardy optimization was success.")
        causality = true
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])

    reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)

    return causality
end
