mutable struct NevanlinnaSolver{T<:Real}
    imags::ImagDomainData{T}          #imaginary domain data
    reals::RealDomainData{T}          #real domain data
    phis::Vector{Complex{T}}          #phis in schur algorithm
    abcd::Array{Complex{T},3}         #continued fractions
    H_max::Int64                      #upper cut off of H
    H_min::Int64                      #lower cut off of H
    H::Int64                          #current value of H
    ab_coeff::Vector{ComplexF64}      #current solution for H
    hardy_matrix::Array{Complex{T},2} #hardy_matrix for H
    iter_tol::Int64                   #upper bound of iteration
    lambda::Float64                   #regularization parameter for second derivative term
    verbose::Bool                       
end

function NevanlinnaSolver(
                  wn          ::Vector{Complex{T}},
                  gw          ::Vector{Complex{T}},
                  N_real      ::Int64,
                  w_max       ::Float64,
                  eta         ::Float64,
                  sum_rule    ::Float64,
                  H_max       ::Int64,
                  iter_tol    ::Int64,
                  lambda      ::Float64
                  ;
                  verbose     ::Bool=false,
                  pick_check  ::Bool=true,
                  optimization::Bool=true,
                  mesh        ::Symbol=:linear,
                  ham_option  ::Bool=false #option for using in Hamburger moment problem
                  )::NevanlinnaSolver{T} where {T<:Real}

    if N_real%2 == 1
        error("N_real must be even number!")
    end

    @assert length(wn) == length(gw)
    N_imag = length(wn) 

    if pick_check
        opt_N_imag =  calc_opt_N_imag(N_imag, wn, gw, verbose=verbose)
    else 
        opt_N_imag = N_imag
    end

    imags = ImagDomainData(wn, gw, opt_N_imag)
    reals = RealDomainData(N_real, w_max, eta, sum_rule, T=T, mesh=mesh)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)

    H_min::Int64 = 1
    ab_coeff = zeros(ComplexF64, 2*H_min)
    hardy_matrix = calc_hardy_matrix(reals, H_min)

    sol = NevanlinnaSolver(imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, verbose)

    if ham_option
        return sol
    end
    
    if optimization
        calc_H_min(sol)
    else
        evaluation!(sol)
    end

    return sol
end

function calc_H_min(sol::NevanlinnaSolver{T},)::Nothing where {T<:Real}
    H_bound::Int64 = 50
    for iH in 1:H_bound
        println("H=$(iH)")
        zero_ab_coeff = zeros(ComplexF64, 2*iH)

        causality, optim = hardy_optim!(sol, iH, zero_ab_coeff, iter_tol=500)

        #break if we find optimal H in which causality is preserved and optimize is successful
        if causality && optim
            sol.H_min = sol.H
            break
        end

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end

        if iH == H_bound
            error("H_min does not exist")
        end
    end
end


function solve!(sol::NevanlinnaSolver{T})::Nothing where {T<:Real}
    ab_coeff  = copy(sol.ab_coeff)

    for iH in sol.H_min:sol.H_max
        println("H=$(iH)")
        causality, optim = hardy_optim!(sol, iH, ab_coeff)

        #break if we face instability of optimization
        if !(causality && optim)
            break
        end

        ab_coeff  = copy(sol.ab_coeff)
        push!(ab_coeff, 0.0+0.0*im)
        push!(ab_coeff, 0.0+0.0*im)

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
end
