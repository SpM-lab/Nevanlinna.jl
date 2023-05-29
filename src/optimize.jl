function calc_functional(
                    sol::NevanlinnaSolver{T},
                    H::Int64, 
                    ab_coeff::Vector{Complex{S}}, 
                    hardy_matrix::Array{Complex{T},2};
                    )::Float64 where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
    green = im * (one(T) .+ theta) ./ (one(T) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = integrate(sol.reals.freq, A)
    second_der = integrate_squared_second_deriv(sol.reals.freq, A) 

    max_theta = findmax(abs.(param))[1]
    func = abs(sol.reals.sum_rule-tot_int)^2 + sol.lambda*second_der

    return func
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
    
    if  !(Optim.converged(res)) && sol.verbose
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

function calc_functional(
                sol         ::HamburgerNevanlinnaSolver{T},
                H           ::Int64, 
                ab_coeff    ::Vector{Complex{S}}, 
                hardy_matrix::Array{Complex{T},2};
                )::Float64 where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    theta = (sol.nev_st.abcd[1,1,:].* param .+ sol.nev_st.abcd[1,2,:]) ./ (sol.nev_st.abcd[2,1,:].*param .+ sol.nev_st.abcd[2,2,:])
    nev_val = im * (one(T) .+ theta) ./ (one(T) .- theta)

    P, Q, G, D = calc_PQGD(sol.mat_real_omega, sol.p, sol.q, sol.gamma, sol.delta)
    val = (- G .- nev_val .* D) ./ (P .+ nev_val .* Q)

    A = Float64.(imag(val)./pi)

    tot_int = integrate(sol.nev_st.reals.freq, A)
    second_der = integrate_squared_second_deriv(sol.nev_st.reals.freq, A) 

    func = abs(sol.nev_st.reals.sum_rule-tot_int)^2 + sol.nev_st.lambda*second_der

    return func
end



function hardy_optim!(
              sol     ::HamburgerNevanlinnaSolver{T},
              H       ::Int64,
              ab_coeff::Array{ComplexF64,1};
              iter_tol::Int64=sol.nev_st.iter_tol,
              )::Tuple{Bool, Bool} where {T<:Real}

    loc_hardy_matrix = calc_hardy_matrix(sol.nev_st.reals, H)

    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(sol, H, x, loc_hardy_matrix)
    end

    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = sol.verbose))
    
    if  !(Optim.converged(res)) && sol.verbose
        println("Faild to optimize!")
    end

    causality = check_causality(loc_hardy_matrix, Optim.minimizer(res), verbose=sol.verbose)

    #causality = evaluation!(reals, abcd, H, Optim.minimizer(res), hardy_matrix, verbose=verbose)
    if causality && (Optim.converged(res))
        sol.nev_st.H = H
        sol.nev_st.ab_coeff = Optim.minimizer(res)
        sol.nev_st.hardy_matrix = loc_hardy_matrix
        hamburger_evaluation!(sol, verbose=false)
    end
    
    return causality, (Optim.converged(res))
end


