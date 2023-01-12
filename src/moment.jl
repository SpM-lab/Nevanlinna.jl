mutable struct HamburgerNevanlinnaSolver{T<:Real}
    moments::Vector{Complex{T}}      
    N_moments_::Int64
    N::Int64
    n1::Int64
    n2::Int64
    isPSD::Bool
    isSingular::Bool
    isProper::Bool
    isDegenerate::Bool
    p::Vector{Complex{T}}      
    q::Vector{Complex{T}}      
    gamma::Vector{Complex{T}}      
    delta::Vector{Complex{T}}      
    hankel::Array{Complex{T},2}
    mat_real_omega::Array{Complex{T},2}
    val::Vector{Complex{T}}      
    nev_struct::NevanlinnaSolver{T}
    verbose::Bool
end

function HamburgerNevanlinnaSolver(
                          moments::Vector{Complex{T}},
                          matsu_omega::Vector{Complex{T}},
                          matsu_green::Vector{Complex{T}},
                          N_real::Int64,
                          omega_max::Float64,
                          eta::Float64,
                          sum_rule::Float64,
                          H_max::Int64,
                          iter_tol::Int64,
                          lambda::Float64
                          ;
                          verbose::Bool=false,
                          pick_check=true,
                          optimization=true,
                          mesh::Symbol=:linear
                          )::HamburgerNevanlinnaSolver{T} where {T<:Real}

    N_moments_ = length(moments)
    if N_moments_ % 2 == 0
        error("invalid moment number. Moment number should be odd.")
    end
    N = div((N_moments_ + 1) , 2)

    #generate hankel matrix
    hankel = Array{Complex{T}}(undef, N, N)
    for i in 1:N, j in 1:N
        hankel[i,j] = moments[i+j-1]
    end

    n1, n2, isDegenerate, isPSD, isSingular, isProper = existence_condition(hankel)

    p, q, gamma, delta = coefficient_lists(moments, hankel, n1, n2, isDegenerate, isPSD, isSingular, isProper)

    if N_real%2 == 1
        error("N_real must be even number!")
    end

    @assert length(matsu_omega) == length(matsu_green)
    N_imag = length(matsu_omega) 

    embed_nev_val = Vector{Complex{T}}(undef, N_imag)
    for i in 1:N_imag
        z::Complex{T} = matsu_omega[i]
        P, Q, G, D = calc_PQGD(z, p, q, gamma, delta)
        nev_val = -matsu_green[i]
        embed_nev_val[i] = (- nev_val * P - G) / (nev_val * Q + D)
    end

    if pick_check
        opt_N_imag =  calc_opt_N_imag(N_imag, matsu_omega, -embed_nev_val, verbose=verbose)
    else 
        opt_N_imag = N_imag
    end

    #generate input data for Schur
    imags = ImagDomainData(matsu_omega, -embed_nev_val, opt_N_imag)
    reals = RealDomainData(N_real, omega_max, eta, sum_rule, T=T, mesh=mesh)
    
    mat_real_omega  = Array{Complex{T}}(undef, N_real, n2+1)
    for i in 1:N_real, j in 1:(n2 + 1)
        mat_real_omega[i,j]  = reals.freq[j]^(i-1)
    end

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    
    val = zeros(Complex{T}, N_real)

    if optimization
        val, reals, H_min, ab_coeff = calc_H_min(p, q, gamma, delta, mat_real_omega, reals, abcd, lambda, verbose)
        hardy_matrix = calc_hardy_matrix(reals, H_min)
    else
        H_min::Int64 = 1
        ab_coeff = zeros(ComplexF64, 2*H_min)
        hardy_matrix = calc_hardy_matrix(reals, H_min)
        evaluation!(reals, abcd, H_min, ab_coeff, hardy_matrix)
        hamburger_evaluation!(p, q, gamma, delta, val, reals, abcd, ab_coeff, hardy_matrix,verbose=verbose)
    end

    nev_sol = NevanlinnaSolver(imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, verbose)
    
    hamburger_sol = HamburgerNevanlinnaSolver(moments, N_moments_, N, n1, n2, isPSD, isSingular, isProper, isDegenerate, p, q, gamma, delta, hankel, mat_real_omega, val, nev_sol, verbose)
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

function solve!(sol::HamburgerNevanlinnaSolver{T})::Nothing where {T<:Real}
    opt_reals = deepcopy(sol.nev_struct.reals)
    ab_coeff  = copy(sol.nev_struct.ab_coeff)
    causality = false
    optim     = false
    
    for iH in sol.nev_struct.H_min:sol.nev_struct.H_max
        println("H=$(iH)")
        hardy_matrix = calc_hardy_matrix(sol.nev_struct.reals, iH)
        opt_val, opt_reals, ab_coeff, causality, optim = Nevanlinna_Schur(sol.p, sol.q, sol.gamma, sol.delta, sol.mat_real_omega, sol.nev_struct.reals, sol.nev_struct.abcd, iH, ab_coeff, hardy_matrix, sol.nev_struct.iter_tol, sol.nev_struct.lambda, sol.nev_struct.verbose)

        #break if we face instability of optimization
        if causality && optim
            sol.val                     = copy(opt_val)
            sol.nev_struct.reals        = deepcopy(opt_reals)
            sol.nev_struct.hardy_matrix = hardy_matrix
            sol.nev_struct.H            = iH
            sol.nev_struct.ab_coeff     = copy(ab_coeff)
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

function existence_condition(hankel::Matrix{Complex{T}}
                            )::Tuple{Int64, Int64, Bool, Bool, Bool, Bool} where {T<:Real}
    N = size(hankel,1)
                
    #compute rank
    n1::Int64 = rank(hankel)
    n2::Int64 = 2*N - n1
    println("Rank of Hankel matrix:$(n1)")

    if n1 == 0
        error("Meeting degenerate 0 matrix.")
    end

    #check degeneracy
    if hankel[1,:] == zeros(Complex{T},N) && hankel[:,1] == zeros(Complex{T},N)
        println("Degenerate")
        isDegenerate = true
    else
    println("Non-degenerate")
    isDegenerate = false
    end

    #check positive semi-definiteness
    PSD_test = hankel .+ T(1e-250).*Matrix{Complex{T}}(I, N, N)
    isPSD = issuccess(cholesky(PSD_test,check = false))
    if isPSD
        println("Postive semi-definite")
    else
        println("Meeting non positive semi-definite matrix in moment calculation.")
    end

    #check singularity
    if n1 < N
        isSingular = true
        println("Singular")
    else
        isSingular = false
        println("Non-singular")
        if isPSD
            println("Positive definite")
        end
    end

    #check properness
    tl_hankel = hankel[1:n1, 1:n1]
    if rank(tl_hankel) < n1 
        isProper = false
        println("Non-proper")
    else
        isProper = true
        println("Proper")
    end

    return n1, n2, isDegenerate, isPSD, isSingular, isProper
end


function coefficient_lists(moments::Vector{Complex{T}},
                           hankel::Matrix{Complex{T}},
                           n1::Int64,
                           n2::Int64,
                           isDegenerate::Bool,
                           isPSD::Bool,
                           isSingular::Bool,
                           isProper::Bool
                          )::Tuple{Vector{Complex{T}}, Vector{Complex{T}}, Vector{Complex{T}}, Vector{Complex{T}}} where {T<:Real}

    N = size(hankel,1)

    #fill extended hankel matrix for calculation
    if !(isSingular)
        extended_hankel = zeros(Complex{T}, N+1, N+1)
        extended_hankel[1:N,1:N] .= hankel
        for i in 1:(N-1)
            extended_hankel[i, N+1] = moments[i+N]
        end
        for j in 1:(N-1)
            extended_hankel[N+1, j] = moments[j+N]
        end
        extended_hankel[N, N+1] = Complex{T}(1.0)
    else
        extended_hankel = copy(hankel)
    end

    #p, q
    p = zeros(Complex{T},n1+1)
    q = zeros(Complex{T},n2+1)
    if isDegenerate
        p[1] = Complex{T}(1.0)
    elseif isProper
        orthogonal_polynomial(extended_hankel, p, n1)
        orthogonal_polynomial(extended_hankel, q, n1 - 1)
    else  #non-proper but not degenerate
        orthogonal_polynomial(extended_hankel, p, n1 - 1)
        #kernel of A_{n2 + 1}
        A = extended_hankel[1:(n1-1),1:(n2+1)]
        q = nullspace(A)[:,1]
        norm::Complex{T} = q[n2]
        for i in 1:n2 
            q[i] /= norm
        end
    end

    #gamma, delta
    sym_p = symmetrizer(p[2:(n1+1)])
    sym_q = symmetrizer(q[2:(n2+1)])
    gamma = sym_p * moments[1:n1]
    delta = sym_q * moments[1:n2]

    return p, q, gamma, delta
end


function symmetrizer(vec::Vector{Complex{T}}
                    )::Matrix{Complex{T}} where {T<:Real}
    dim::Int64 = length(vec)
    mat = Array{Complex{T}}(undef,dim,dim)
    for i in 1:dim
        for j in 1:(dim-i+1)
            mat[i,j] = vec[i+j-1]
        end
        for j in (dim-i+2):dim
            mat[i,j] = Complex{T}(0.0)
        end
    end
    return mat
end

function removeColumn(matrix::Matrix{Complex{T}}, colToRemove::Int64
                     )::Matrix{Complex{T}} where {T<:Real}

    numCols::Int64 = size(matrix,2)

    if colToRemove == 1
        rem_matrix = matrix[:,2:numCols]
    elseif colToRemove == numCols
        rem_matrix = matrix[:,1:numCols-1]
    else
        leftmat    = matrix[:,1:colToRemove-1]
        rightmat   = matrix[:,colToRemove+1:numCols]
        rem_matrix = hcat(leftmat,rightmat)
    end

    return rem_matrix
end

function orthogonal_polynomial(mat::Matrix{Complex{T}},
                               vec::Vector{Complex{T}},
                               order::Int64
                              )::Nothing where {T<:Real}

    sliced_mat = mat[1:order, 1:(order+1)]
    #get cofactor of sliced matrix, as coefficients of the polynomial vector
    for i in 1:(order+1)
        temp = copy(sliced_mat)
        temp = removeColumn(temp, i)
        vec[i] = (-1)^(i+order) * det(temp)
    end
    norm::Complex{T} = vec[order+1]
    for i in 1:(order+1)
        vec[i] /= norm
    end
end

function calc_PQGD(z::Complex{T},
                   p::Vector{Complex{T}},
                   q::Vector{Complex{T}},
                   gamma::Vector{Complex{T}},
                   delta::Vector{Complex{T}}
                   )::Tuple{Complex{T},Complex{T},Complex{T},Complex{T}} where {T<:Real}

    n2::Int64 = length(delta)
    poly_val = Vector{Complex{T}}(undef, (n2 + 1))

    for i in 1:(n2 + 1)
        poly_val[i] = z^(i-1)
    end

    P = sum(poly_val[1:length(p)] .* p)
    Q = sum(poly_val[1:length(q)] .* q)
    G = sum(poly_val[1:length(gamma)] .* gamma)
    D = sum(poly_val[1:length(delta)] .* delta)

    return P, Q, G, D
end

function calc_PQGD(matz::Matrix{Complex{T}},
                   p::Vector{Complex{T}},
                   q::Vector{Complex{T}},
                   gamma::Vector{Complex{T}},
                   delta::Vector{Complex{T}}
                   )::Tuple{Vector{Complex{T}},Vector{Complex{T}},Vector{Complex{T}},Vector{Complex{T}}} where {T<:Real}

    n2::Int64 = length(delta)

    P = matz[:,1:length(p)] * p
    Q = matz[:,1:length(q)] * q
    G = matz[:,1:length(gamma)] * gamma
    D = matz[:,1:length(delta)] * delta

    return P, Q, G, D
end
