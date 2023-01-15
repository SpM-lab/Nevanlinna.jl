function calc_opt_N_imag(N::Int64,
                         matsu_omega::Array{Complex{T},1},
                         matsu_green::Array{Complex{T},1};
                         verbose::Bool=false
                         )::Int64 where {T<:Real}
    @assert N == length(matsu_omega)
    @assert N == length(matsu_green)

    freq = (matsu_omega  .- im) ./ (matsu_omega  .+ im)
    val  = (-matsu_green .- im) ./ (-matsu_green .+ im)

    k::Int64 = 0
    success::Bool = true

    while success
        k += 1
        Pick = Array{Complex{T}}(undef, k, k)

        for j in 1:k
            for i in 1:k
                num = one(T) - val[i]  * conj(val[j])
                den = one(T) - freq[i] * conj(freq[j])
                Pick[i,j] = num / den
            end

            Pick[j,j] += T(1e-250)
        end

        success = issuccess(cholesky(Pick,check = false))

        if k == N
            break
        end
    end

    if verbose
        println("N_imag is setted as $(k-1)")
    end

    return (k-1)
end

function hardy_optim!(
                sol::NevanlinnaSolver{T},
                H::Int64,
                ab_coeff::Array{ComplexF64,1};
                iter_tol::Int64=sol.iter_tol,
                )::Tuple{Bool, Bool} where {T<:Real}

    loc_hardy_matrix = calc_hardy_matrix(sol.reals, H)

    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(sol, H, x, loc_hardy_matrix)
    end

    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = sol.verbose))
    
    if  !(Optim.converged(res))
        println("Faild to optimize!")
    end
    
    causality = check_causality(loc_hardy_matrix, Optim.minimizer(res), verbose=sol.verbose)

    if causality && (Optim.converged(res))
        sol.H = H
        sol.ab_coeff = Optim.minimizer(res)
        sol.hardy_matrix = loc_hardy_matrix
        evaluation!(sol, verbose=false)
    end
    
    return causality, (Optim.converged(res))
end
