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


