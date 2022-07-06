function calc_opt_N_imag(N::Int64,
                         matsu_omega::Array{Complex{T},1},
                         matsu_green::Array{Complex{T},1}
                         )::Int64 where {T<:Real}
    @assert N == length(matsu_omega)
    @assert N == length(matsu_green)

    freq = (matsu_omega .- im) ./ (matsu_omega .+ im)
    val  = (-matsu_green .- im) ./ (-matsu_green .+ im)

    k::Int64 = 0
    success::Bool = true

    while success
        k += 1
        Pick = Array{Complex{T}}(undef, k, k)

        for j in 1:k
            for i in 1:k
                nom = one(T) - val[i] * conj(val[j])
                den = one(T) - freq[i] * conj(freq[j])
                Pick[i,j] = nom / den
            end
            Pick[j,j] += T(1e-250)
        end

        success = issuccess(cholesky(Pick,check = false))

        if k == N
            break
        end
    end

    println("N_imag is setted as $(k-1)")

    return (k-1)
end

 
function Nevanlinna_Schur(reals::RealDomainData{T},
                          abcd::Array{Complex{T},3},
                          H::Int64,
                          ab_coeff::Array{ComplexF64,1},
                          hardy_matrix::Array{Complex{T},2},
                          iter_tol::Int64,
                          lambda::Float64,
                          verbose::Bool=false
                          )::Tuple{RealDomainData{T}, Array{ComplexF64,1}, Bool, Bool} where {T<:Real}

    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(reals, abcd, H, x, hardy_matrix, lambda=lambda)
    end
    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = verbose))
    
    if  !(Optim.converged(res))
        #error("Faild to optimize!")
        println("Faild to optimize!")
    end
    
    causality = evaluation!(reals, abcd, H, Optim.minimizer(res), hardy_matrix)
    
    return reals, Optim.minimizer(res), causality, (Optim.converged(res))
end

function calc_error(reals::RealDomainData{T},
                    matsu_omega::Array{Complex{T},1},
                    matsu_green::Array{Complex{T},1}) where {T<:Real}
        delta_green = Array{Complex{T}}(undef, length(matsu_omega))

        for i in 1:length(matsu_omega)
            int_green = imag.(reals.val)/pi ./ (matsu_omega[i] .- (reals.freq))
            back_i = integrate(real.(reals.freq), int_green)
            delta_green[i] = back_i - matsu_green[i]
        end

        return findmax(abs.(delta_green))[1]
end


function calc_error(imags::ImagDomainData{T},
                    matsu_omega::Vector{Complex{T}},
                    matsu_green::Vector{Complex{T}},
                    phis::Vector{Complex{T}},
                    H::Int64,
                    ab_coeff::Vector{ComplexF64}
                    ) where {T<:Real}

    matsu_abcd = Array{Complex{T}}(undef, 2, 2,  length(matsu_omega))
    for i in 1:length(matsu_omega)
        result = Matrix{Complex{T}}(I, 2, 2)
        z = matsu_omega[i]
        for j in 1:imags.N_imag
            prod = Array{Complex{T}}(undef, 2, 2)
            prod[1,1] = (z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[2,2] = one(T)
            result *= prod
        end
        matsu_abcd[:,:,i] .= result
    end


    matsu_hardy_matrix = Array{Complex{T}}(undef, length(matsu_omega), 2*H)
    for k in 1:H
        matsu_hardy_matrix[:,2*k-1] .=      hardy_basis.(matsu_omega,k-1)
        matsu_hardy_matrix[:,2*k]   .= conj(hardy_basis.(matsu_omega,k-1))
    end
    
    matsu_param = matsu_hardy_matrix*ab_coeff

    matsu_theta = (matsu_abcd[1,1,:].* matsu_param .+ matsu_abcd[1,2,:]) ./ (matsu_abcd[2,1,:].*matsu_param .+ matsu_abcd[2,2,:])

    back_matsu_green = -im * (one(T) .+ matsu_theta) ./ (one(T) .- matsu_theta)

    return findmax(abs.(back_matsu_green - matsu_green))[1]
end

    
