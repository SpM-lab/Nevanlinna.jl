function calc_H_min(p::Vector{Complex{T}},
                    q::Vector{Complex{T}},
                    gamma::Vector{Complex{T}},
                    delta::Vector{Complex{T}},
                    mat_real_omega::Array{Complex{T},2},
                    reals::RealDomainData,
                    abcd::Array{Complex{T},3},
                    lambda::Float64,
                    verbose::Bool=false
                    )::Tuple{Vector{Complex{T}}, RealDomainData{T}, Int64, Array{ComplexF64,1}} where {T<:Real}
    for iH in 1:50
        println("H=$(iH)")
        zero_ab_coeff = zeros(ComplexF64, 2*iH)
        hardy_matrix = calc_hardy_matrix(reals, iH)
        opt_val, opt_reals, ab_coeff, causality, optim = Nevanlinna_Schur(p, q, gamma, delta, mat_real_omega, reals, abcd, iH, zero_ab_coeff, hardy_matrix, 500, lambda, verbose)

        #break if we find optimal H in which causality is preserved and optimize is successful
        if causality && optim
            return opt_val, opt_reals, iH, ab_coeff
            break
        end

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
    error("H_min does not exist")
end


function Nevanlinna_Schur(p::Vector{Complex{T}},
                          q::Vector{Complex{T}},
                          gamma::Vector{Complex{T}},
                          delta::Vector{Complex{T}},
                          mat_real_omega::Array{Complex{T},2},
                          reals::RealDomainData{T},
                          abcd::Array{Complex{T},3},
                          H::Int64,
                          ab_coeff::Array{ComplexF64,1},
                          hardy_matrix::Array{Complex{T},2},
                          iter_tol::Int64,
                          lambda::Float64,
                          verbose::Bool=false
                          )::Tuple{Vector{Complex{T}}, RealDomainData{T}, Array{ComplexF64,1}, Bool, Bool} where {T<:Real}

    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(p, q, gamma, delta, mat_real_omega, reals, abcd, H, x, hardy_matrix, lambda=lambda)
    end

    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = verbose))
    
    if  !(Optim.converged(res))
        println("Faild to optimize!")
    end
    
    causality = evaluation!(reals, abcd, H, Optim.minimizer(res), hardy_matrix, verbose=verbose)
    val = zeros(Complex{T}, reals.N_real)
    hamburger_evaluation!(p, q, gamma, delta, val, reals, abcd, Optim.minimizer(res), hardy_matrix,verbose=verbose)
    
    return val, reals, Optim.minimizer(res), causality, (Optim.converged(res))
end

function calc_functional(p::Vector{Complex{T}},
                         q::Vector{Complex{T}},
                         gamma::Vector{Complex{T}},
                         delta::Vector{Complex{T}},
                         mat_real_omega::Array{Complex{T},2},
                         reals::RealDomainData{T}, 
                         abcd::Array{Complex{T},3}, 
                         H::Int64, 
                         ab_coeff::Vector{Complex{S}}, 
                         hardy_matrix::Array{Complex{T},2};
                         lambda::Float64 = 1e-5
                         )::Float64 where {S<:Real, T<:Real}
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    nev_val = im * (one(T) .+ theta) ./ (one(T) .- theta)

    val = zeros(Complex{T}, reals.N_real)
#=
    for i in 1:reals.N_real
        z::Complex{T} = reals.freq[i]
        P, Q, G, D = calc_PQGD(z, p, q, gamma, delta)
        val[i] = (- G - nev_val[i] * D) / (P + nev_val[i] * Q)
    end
=#
    P, Q, G, D = calc_PQGD(mat_real_omega, p, q, gamma, delta)
    val = (- G .- nev_val .* D) ./ (P .+ nev_val .* Q)

    A = Float64.(imag(val)./pi)

    tot_int = integrate(reals.freq, A)
    second_der = integrate_squared_second_deriv(reals.freq, A) 

    max_theta = findmax(abs.(param))[1]
    func = abs(reals.sum-tot_int)^2 + lambda*second_der

    return func
end

function hamburger_evaluation!(p::Vector{Complex{T}},
                               q::Vector{Complex{T}},
                               gamma::Vector{Complex{T}},
                               delta::Vector{Complex{T}},
                               val::Vector{Complex{T}},
                               reals::RealDomainData{T}, 
                               abcd::Array{Complex{T},3}, 
                               ab_coeff::Vector{Complex{S}}, 
                               hardy_matrix::Array{Complex{T},2};
                               verbose::Bool=false
                               )::Bool where {S<:Real, T<:Real}

    param = hardy_matrix * ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        #if verbose
        #    println("max_theta=",max_theta)
        #    println("hardy optimization was success.")
        #end
        causality = true

        theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])

        reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)

        for i in 1:reals.N_real
            z::Complex{T} = reals.freq[i]
            P, Q, G, D = calc_PQGD(z, p, q, gamma, delta)
            val[i] = (- G - reals.val[i] * D) / (P + reals.val[i] * Q)
        end
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end

    return causality
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
                    H::Int64, 
                    ab_coeff::Vector{Complex{S}}, 
                    hardy_matrix::Array{Complex{T},2};
                    verbose::Bool=false
                    )::Bool where {S<:Real, T<:Real}

    causality = check_causality(hardy_matrix, ab_coeff, verbose=verbose)
    param = hardy_matrix*ab_coeff
    if causality
        theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
        reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)
    end
#=
    param = hardy_matrix*ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        if verbose
           println("max_theta=",max_theta)
           println("hardy optimization was success.")
        end
        causality = true

        theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
        reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end
=#

    return causality
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
        println("Faild to optimize!")
    end
    
    causality = evaluation!(reals, abcd, H, Optim.minimizer(res), hardy_matrix, verbose=verbose)
    
    return reals, Optim.minimizer(res), causality, (Optim.converged(res))
end

function calc_H_min(reals::RealDomainData,
                    abcd::Array{Complex{T},3},
                    lambda::Float64,
                    verbose::Bool=false
                    )::Tuple{RealDomainData{T}, Int64, Array{ComplexF64,1}} where {T<:Real}
    for iH in 1:50
        println("H=$(iH)")
        zero_ab_coeff = zeros(ComplexF64, 2*iH)
        hardy_matrix = calc_hardy_matrix(reals, iH)
        opt_reals, ab_coeff, causality, optim = Nevanlinna_Schur(reals, abcd, iH, zero_ab_coeff, hardy_matrix, 500, lambda, verbose)

        #break if we find optimal H in which causality is preserved and optimize is successful
        if causality && optim
            return opt_reals, iH, ab_coeff
            break
        end

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
    error("H_min does not exist")
end


