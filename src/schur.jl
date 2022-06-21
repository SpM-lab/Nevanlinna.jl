function Nevanlinna_Schur(N_imag::Int64, 
                    omega::Array{T,1}, 
                    green::Array{Complex{T},1},
                    N_real::Int64,
                    omega_max::Float64,
                    eta::Float64,
                    H::Int64,
                    verbose::Bool=false)::Tuple{ImagDomainData{T}, RealDomainData{T}} where {T<:Real}
    if N_real%2 == 1
        error("N_real must be even number!")
    end
    
    imags = ImagDomainData(N_imag, omega, green)
    reals     = RealDomainData(N_real, omega_max, eta, T)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    hardy_matrix = calc_hardy_matrix(reals, H)
    
    ab_coeff  = zeros(Complex{T}, 2*H) 
    
    functional = x->Nevanlinna.calc_functional(reals, abcd, H, x, hardy_matrix)
    
    function jacobian(J::Vector, x)
        @assert length(J) == length(x)
        J .= gradient(functional, x)[1] 
    end
    J = similar(ab_coeff)
   
    if verbose
        res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                        Optim.Options(iterations = 100000,
                                      show_trace = true))
    else 
        res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                        Optim.Options(iterations = 100000,
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


function calc_opt_N_imag(N::Int64,
                         omega::Array{T,1},
                         green::Array{Complex{T},1})::Int64 where {T<:Real}
    @assert N == length(omega)
    @assert N == length(green)

    freq = (omega*im .- im) ./ (omega*im .+ im)
    val  = (-green .- im) ./ (-green .+ im)

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
