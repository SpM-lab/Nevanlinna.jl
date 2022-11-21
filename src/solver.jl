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

function NevanlinnaSolver(N_imag::Int64,
                          matsu_omega::Vector{Complex{T}},
                          matsu_green::Vector{Complex{T}},
                          N_real::Int64,
                          omega_max::Float64,
                          eta::Float64,
                          sum::Float64,
                          H_max::Int64,
                          iter_tol::Int64,
                          lambda::Float64
                          ;
                          verbose::Bool=false,
                          pick_check=true,
                          optimization=true,
                          mesh::Symbol=:linear
                          )::NevanlinnaSolver{T} where {T<:Real}
    if N_real%2 == 1
        error("N_real must be even number!")
    end

    if pick_check
        opt_N_imag =  calc_opt_N_imag(N_imag, matsu_omega, matsu_green, verbose=verbose)
    else 
        opt_N_imag = N_imag
    end

    imags = ImagDomainData(matsu_omega, matsu_green, opt_N_imag)
    reals = RealDomainData(N_real, omega_max, eta, sum, T=T, mesh=mesh)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    
    if optimization
        reals, H_min, ab_coeff = calc_H_min(reals, abcd, lambda, verbose)
        hardy_matrix = calc_hardy_matrix(reals, H_min)
    else
        H_min::Int64 = 1
        ab_coeff = zeros(ComplexF64, 2*H_min)
        hardy_matrix = calc_hardy_matrix(reals, H_min)
        evaluation!(reals, abcd, H_min, ab_coeff, hardy_matrix)
    end

    sol = NevanlinnaSolver(imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, verbose)
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

function solve!(sol::NevanlinnaSolver{T})::Nothing where {T<:Real}
    opt_reals = deepcopy(sol.reals)
    ab_coeff  = copy(sol.ab_coeff)
    causality = false
    optim     = false
    
    for iH in sol.H_min:sol.H_max
        println("H=$(iH)")
        hardy_matrix = calc_hardy_matrix(sol.reals, iH)
        opt_reals, ab_coeff, causality, optim = Nevanlinna_Schur(sol.reals, sol.abcd, iH, ab_coeff, hardy_matrix, sol.iter_tol, sol.lambda, sol.verbose)

        #break if we face instability of optimization
        if causality && optim
            sol.reals        = deepcopy(opt_reals)
            sol.hardy_matrix = hardy_matrix
            sol.H            = iH
            sol.ab_coeff     = copy(ab_coeff)
        else
            break
        end

        push!(ab_coeff, 0.0+0.0*im)
        push!(ab_coeff, 0.0+0.0*im)

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
end
