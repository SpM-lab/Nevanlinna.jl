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

function check_causality(hardy_matrix::Array{Complex{T},2},
                         ab_coeff::Vector{Complex{S}};
                         verbose::Bool=false
                         )::Bool where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        if verbose
           println("max_theta=",max_theta)
           println("hardy optimization was success.")
        end
        causality = true
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end
    return causality
end

function evaluation!(sol::NevanlinnaSolver{T};
                     verbose::Bool=false
                    )::Bool where {T<:Real}

    causality = check_causality(sol.hardy_matrix, sol.ab_coeff, verbose=verbose)
    if causality
        param = sol.hardy_matrix*sol.ab_coeff
        theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
        sol.reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)
    end

    return causality
end

function hamburger_evaluation!(
                sol    ::HamburgerNevanlinnaSolver{T};
                verbose::Bool=false
                )::Bool where {T<:Real}

    causality = check_causality(sol.nev_st.hardy_matrix, sol.nev_st.ab_coeff, verbose=verbose)

    if causality
        param = sol.nev_st.hardy_matrix*sol.nev_st.ab_coeff
        theta = (sol.nev_st.abcd[1,1,:].* param .+ sol.nev_st.abcd[1,2,:]) ./ (sol.nev_st.abcd[2,1,:].*param .+ sol.nev_st.abcd[2,2,:])

        sol.nev_st.reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)

        P, Q, G, D = calc_PQGD(sol.mat_real_omega, sol.p, sol.q, sol.gamma, sol.delta)
        sol.val .= (- G .- sol.nev_st.reals.val .* D) ./ (P .+ sol.nev_st.reals.val .* Q)
    end

    return causality
end


