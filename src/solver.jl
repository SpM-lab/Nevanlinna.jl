mutable struct NevanlinnaSolver{T<:Real}
    matsu_omega::Vector{Complex{T}}
    matsu_green::Vector{Complex{T}}
    imags::ImagDomainData{T}
    reals::RealDomainData{T}
    phis::Vector{Complex{T}}
    abcd::Array{Complex{T},3}
    H_max::Int64
    H_min::Int64
    H::Int64
    ab_coeff::Vector{ComplexF64}
    hardy_matrix::Array{Complex{T},2}
    lambda::Float64
    verbose::Bool
end

function NevanlinnaSolver(N_imag::Int64,
                          matsu_omega::Vector{Complex{T}},
                          matsu_green::Vector{Complex{T}},
                          N_real::Int64,
                          omega_max::Float64,
                          eta::Float64,
                          H_max::Int64,
                          lambda::Float64,
                          verbose::Bool=false)::NevanlinnaSolver{T} where {T<:Real}
    if N_real%2 == 1
        error("N_real must be even number!")
    end

    imags = ImagDomainData(matsu_omega, matsu_green, N_imag)
    reals = RealDomainData(N_real, omega_max, eta, T=T)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    
    H_min = calc_H_min(reals, abcd, lambda, true)

    ab_coeff  = zeros(ComplexF64, 2*H_min)
    hardy_matrix = calc_hardy_matrix(reals, H_min)

    evaluation(reals, abcd, H_min, ab_coeff, hardy_matrix)

    sol = NevanlinnaSolver(matsu_omega, matsu_green, imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, lambda, verbose)
end

function calc_H_min(reals::RealDomainData,
                    abcd::Array{Complex{T},3},
                    lambda::Float64,
                    verbose::Bool=false
                    )::Int64 where {T<:Real}
    for iH in 1:50
        println("H= $(iH)")
        ab_coeff = zeros(ComplexF64, 2*iH)
        hardy_matrix = calc_hardy_matrix(reals, iH)
        result = Nevanlinna_Schur(reals, abcd, iH, ab_coeff, hardy_matrix, 500, lambda, verbose)
        causality = result[3]
        optim = result[4]

        if causality && optim
            return iH
            break
        end

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
end

function solve(sol::NevanlinnaSolver) where {T<:Real}
    pre_error = 1.0
    error     = 0.0
    ab_coeff  = zeros(ComplexF64, 2*sol.H_min)
    for iH in sol.H_min:sol.H_max
        println("H= $(iH)")
        hardy_matrix = calc_hardy_matrix(sol.reals, iH)
        result = Nevanlinna_Schur(sol.reals, sol.abcd, iH, ab_coeff, hardy_matrix, 100000, sol.lambda, sol.verbose)
        opt_real  = result[1]
        ab_coeff .= result[2]

        error = calc_error(opt_real, sol.matsu_omega, sol.matsu_green)

        if error > pre_error
            break
        else
            sol.reals        = opt_real
            pre_error        = error
            sol.hardy_matrix = hardy_matrix
            sol.H            = iH
            sol.ab_coeff     = copy(ab_coeff)
        end

        push!(ab_coeff, 0.0+0.0*im)
        push!(ab_coeff, 0.0+0.0*im)

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
end
