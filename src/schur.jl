<<<<<<< HEAD
=======
function Nevanlinna_Schur(N_imag::Int64, 
                    omega::Array{T,1}, 
                    green::Array{Complex{T},1},
                    N_real::Int64,
                    omega_max::Float64,
                    eta::Float64,
                    H::Int64;
                    verbose::Bool=false,
                    iterations = 100000,
                    lambda = 1e-5
                    )::Tuple{ImagDomainData{T}, RealDomainData{T}} where {T<:Real}
    if N_real%2 == 1
        error("N_real must be even number!")
    end
    
    imags = ImagDomainData(N_imag, omega, green)
    reals = RealDomainData(N_real, omega_max, eta, T=T)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    hardy_matrix = calc_hardy_matrix(reals, H)
    
    ab_coeff  = zeros(ComplexF64, 2*H) 
    
    function functional(x::Vector{ComplexF64})::Float64
        return Nevanlinna.calc_functional(reals, abcd, H, x, hardy_matrix, lambda=lambda)
    end
    
    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end
   
    if verbose
        res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                        Optim.Options(iterations = iterations,
                                      show_trace = true))
    else 
        res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                        Optim.Options(iterations = iterations,
                                      show_trace = false))
    end

    #=
    if verbose
        res = optimize(functional, jacobian, ab_coeff, ConjugateGradient(), 
                        Optim.Options(iterations = 100000,
                                      show_trace = true))
    else 
        res = optimize(functional, jacobian, ab_coeff, ConjugateGradient(), 
                        Optim.Options(iterations = 100000,
                                      show_trace = false))
    end
    =#
 
    
    if  !(Optim.converged(res))
        error("Faild to optimize!")
    end
    
    evaluation(reals, abcd, H, Optim.minimizer(res), hardy_matrix)
    
    return imags, reals
end


>>>>>>> 974716fd2b792a13bb2b48f2197ac7b1aad41246
function calc_opt_N_imag(N::Int64,
                         matsu_omega::Array{Complex{T},1},
                         matsu_green::Array{Complex{T},1})::Int64 where {T<:Real}
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
                          verbose::Bool=false)::Tuple{RealDomainData{T}, Array{ComplexF64,1}, Bool, Bool} where {T<:Real}

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
    
    causality = evaluation(reals, abcd, H, Optim.minimizer(res), hardy_matrix)
    
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
