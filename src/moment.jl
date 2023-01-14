mutable struct HamburgerNevanlinnaSolver{T<:Real}
    moments       ::Vector{Complex{T}}      
    N_moments_    ::Int64
    N             ::Int64
    n1            ::Int64
    n2            ::Int64
    isPSD         ::Bool
    isSingular    ::Bool
    isProper      ::Bool
    isDegenerate  ::Bool
    p             ::Vector{Complex{T}}      
    q             ::Vector{Complex{T}}      
    gamma         ::Vector{Complex{T}}      
    delta         ::Vector{Complex{T}}      
    hankel        ::Array{Complex{T},2}
    mat_real_omega::Array{Complex{T},2}
    val           ::Vector{Complex{T}}      
    nev_st        ::NevanlinnaSolver{T}
    verbose       ::Bool
end

function HamburgerNevanlinnaSolver(
                    moments     ::Vector{Complex{T}},
                    matsu_omega ::Vector{Complex{T}},
                    matsu_green ::Vector{Complex{T}},
                    N_real      ::Int64,
                    omega_max   ::Float64,
                    eta         ::Float64,
                    sum_rule    ::Float64,
                    H_max       ::Int64,
                    iter_tol    ::Int64,
                    lambda      ::Float64
                    ;
                    verbose     ::Bool=false,
                    pick_check  ::Bool=true,
                    optimization::Bool=true,
                    mesh        ::Symbol=:linear
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

    H_min::Int64 = 1
    ab_coeff = zeros(ComplexF64, 2*H_min)
    hardy_matrix = calc_hardy_matrix(reals, H_min)

    nev_sol = NevanlinnaSolver(imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, verbose)
    ham_nev_sol = HamburgerNevanlinnaSolver(moments, N_moments_, N, n1, n2, isPSD, isSingular, isProper, isDegenerate, p, q, gamma, delta, hankel, mat_real_omega, val, nev_sol, verbose)

    if optimization
        calc_H_min(ham_nev_sol)
    else
        hamburger_evaluation!(ham_nev_sol)
    end

    return ham_nev_sol
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

    max_theta = findmax(abs.(param))[1]
    func = abs(sol.nev_st.reals.sum-tot_int)^2 + sol.nev_st.lambda*second_der

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
        #return calc_functional(sol.p, sol.q, sol.gamma, sol.delta, sol.mat_real_omega, sol.nev_st.reals, sol.nev_st.abcd, H, x, loc_hardy_matrix, lambda=sol.nev_st.lambda)
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

    #causality = evaluation!(reals, abcd, H, Optim.minimizer(res), hardy_matrix, verbose=verbose)
    if causality && (Optim.converged(res))
        sol.nev_st.H = H
        sol.nev_st.ab_coeff = Optim.minimizer(res)
        sol.nev_st.hardy_matrix = loc_hardy_matrix
#        evaluation!(sol.nev_st, verbose=false)
        hamburger_evaluation!(sol, verbose=false)
    end
    
    return causality, (Optim.converged(res))
end


function calc_H_min(sol::HamburgerNevanlinnaSolver{T})::Nothing where {T<:Real}
    H_bound::Int64 = 50
    for iH in 1:H_bound
        println("H=$(iH)")
        zero_ab_coeff = zeros(ComplexF64, 2*iH)

        causality, optim = hardy_optim!(sol, iH, zero_ab_coeff, iter_tol=500)

        #break if we find optimal H in which causality is preserved and optimize is successful
        if causality && optim
            sol.nev_st.H_min = sol.nev_st.H
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


function solve!(sol::HamburgerNevanlinnaSolver{T})::Nothing where {T<:Real}
    ab_coeff  = copy(sol.nev_st.ab_coeff)
    
    for iH in sol.nev_st.H_min:sol.nev_st.H_max
        println("H=$(iH)")
        causality, optim = hardy_optim!(sol, iH, ab_coeff)

        #break if we face instability of optimization
        if !(causality && optim)
            break
        end

        ab_coeff = copy(sol.nev_st.ab_coeff)
        push!(ab_coeff, 0.0+0.0*im)
        push!(ab_coeff, 0.0+0.0*im)

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end
    end
end

function hamburger_evaluation!(
                sol    ::HamburgerNevanlinnaSolver{T};
                verbose::Bool=false
                )::Bool where {T<:Real}

    causality = check_causality(sol.nev_st.hardy_matrix, sol.nev_st.ab_coeff, verbose=verbose)

    if causality
        param = sol.nev_st.hardy_matrix*sol.nev_st.ab_coeff
        theta = (sol.nev_st.abcd[1,1,:].* param .+ sol.nev_st.abcd[1,2,:]) ./ (sol.nev_st.abcd[2,1,:].*param .+ sol.nev_st.abcd[2,2,:])

        sol.nev_st.reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)

        P, Q, G, D = calc_PQGD(sol.mat_real_omega, sol.p, sol.q, sol.gamma, sol.delta)
        sol.val .= (- G .- sol.nev_st.reals.val .* D) ./ (P .+ sol.nev_st.reals.val .* Q)
    else
        println("hardy optimization was failure.")
    end

    return causality
end

function existence_condition(
                hankel::Matrix{Complex{T}}
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


function coefficient_lists(
                  moments     ::Vector{Complex{T}},
                  hankel      ::Matrix{Complex{T}},
                  n1          ::Int64,
                  n2          ::Int64,
                  isDegenerate::Bool,
                  isPSD       ::Bool,
                  isSingular  ::Bool,
                  isProper    ::Bool
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

function removeColumn(
                    matrix     ::Matrix{Complex{T}}, 
                    colToRemove::Int64
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

function orthogonal_polynomial(
                    mat  ::Matrix{Complex{T}},
                    vec  ::Vector{Complex{T}},
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

function calc_PQGD(z    ::Complex{T},
                   p    ::Vector{Complex{T}},
                   q    ::Vector{Complex{T}},
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

function calc_PQGD(matz ::Matrix{Complex{T}},
                   p    ::Vector{Complex{T}},
                   q    ::Vector{Complex{T}},
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

