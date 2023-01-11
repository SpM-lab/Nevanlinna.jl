mutable struct Hamburger_NevanlinnaSolver{T<:Real}
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
    nev_struct::NevanlinnaSolver{T}
end

function Hamburger_NevanlinnaSolver(
                          moments::Vector{Complex{T}},
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
                          )::Hamburger_NevanlinnaSolver{T} where {T<:Real}
    N_moments_ = length(moments)
    if N_moments%2 == 0
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


#kokomade
    if N_real%2 == 1
        error("N_real must be even number!")
    end

    @assert length(matsu_omega) == length(matsu_green)
    N_imag = length(matsu_omega) 

    if pick_check
        opt_N_imag =  calc_opt_N_imag(N_imag, matsu_omega, matsu_green, verbose=verbose)
    else 
        opt_N_imag = N_imag
    end

    #generate input data for Schur
    orig_imags = ImagDomainData(matsu_omega, matsu_green, opt_N_imag)
    embed_imags = deepcopy(orig_imags)
    reals = RealDomainData(N_real, omega_max, eta, sum, T=T, mesh=mesh)

    poly_val_imag = Vector{Complex{T}}(undef, (n2 + 1))

    for i in 1:orig_imags.N_imag
      z::Complex{T} = orig_imags.freq[i]
      for j in 1:(n2 + 1)
        poly_val_imag[j] = pow(z, j-1)
      end
      G = poly_val_imag[1:length(gamma)] * gamma
      D = poly_val_imag[1:length(delta)] * delta
      P = poly_val_imag[1:length(p)] * p
      Q = poly_val_imag[1:length(q)] * q
      nev_val = im * (one(T) + orig_imags.val[i]) / (one(T) - orig_imags.val[i])
      embed_imags.val[i] = ( (- nev_val * P - G) / (nev_val * Q + D) )
    end

    phis = calc_phis(embed_imags)
    abcd = calc_abcd(embed_imags, reals, phis)
    
    if optimization
        reals, H_min, ab_coeff = calc_H_min(reals, abcd, lambda, verbose)
        hardy_matrix = calc_hardy_matrix(reals, H_min)
    else
        H_min::Int64 = 1
        ab_coeff = zeros(ComplexF64, 2*H_min)
        hardy_matrix = calc_hardy_matrix(reals, H_min)
        evaluation!(reals, abcd, H_min, ab_coeff, hardy_matrix)
    end

    nev_sol = NevanlinnaSolver(embed_imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, verbose)
    
    hamburger_sol = Hamburger_NevanlinnaSolver(moments, N_moments_, N, n1, n2, isPSD, isSingular, isProer, isDegenerate, p, q, gamma, elta, hankel, nev_sol)
end

function solve!(sol::Hamburger_NevanlinnaSolver{T})::Nothing where {T<:Real}

    param = sol.nev_struct.hardy_matrix*sol.nev_struct.ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        if verbose
           println("max_theta=",max_theta)
           println("hardy optimization was success.")
        end
        causality = true

        theta = (sol.nev_struct.abcd[1,1,:].* param .+ sol.nev_struct.abcd[1,2,:]) ./ (sol.nev_struct.abcd[2,1,:].*param .+ sol.nev_struct.abcd[2,2,:])
        sol.nev_struct.reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end
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
    isDegenerate::Bool = true
  else
    println("Non-degenerate")
    isDegenerate::Bool = false
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
    isSingular::Bool = true
    println("Singular")
  else
    isSingular::Bool = false
    println("Non-singular")
    if isPSD
      println("Positive definite")
    end
  end

  #check properness
  tl_hankel = hankel[1:n1, 1:n1]
  if rank(tl_hankel) < n1 
    isProper::Bool = false
    println("Non-proper")
  else
    isProper::Bool = true
    println("Proper")
  end

  return (n1,n2,isDegenerate,isPSD,isSingular,isProper)
end


function coefficient_lists(moments::Vector{Complex{T}},
                           hankel::Matrix{Complex{T}},
                           n1::Int64,
                           n2::Int64,
                           isDegenerate::Bool,
                           isPSD::Bool,
                           isSingular::Bool,
                           isProper::Bool
                          )::Tuple{Vector{Complex{T}}, {Vector{Complex{T}}, {Vector{Complex{T}}, {Vector{Complex{T}} where {T<:Real}

  N = size(hankel,1)

  #fill extended hankel matrix for calculation
  if ! isSingular
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
  else if isProper
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
  gamma = zeros(Complex{T},n1)
  delta = zeros(Complex{T},n2)
  sym_p = symmetrizer(p[2:n1])
  sym_q = symmetrizer(q[2:n2])
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
  else if colToRemove == numCols
    rem_matrix = matrix[:,1:numCols-1]
  else
    leftmat    = matrix[:,1:colToRemove-1]
    rightmat   = matrix[:,colToRemove+1:numCols]
    rem_matrix = hcat(leftmat,rightmat)
  end

  return rem_matrix
}

functionvoid orthogonal_polynomial(mat::Matrix{Complex{T}},
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
