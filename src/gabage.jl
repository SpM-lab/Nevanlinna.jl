function calc_functional(reals::RealDomainData)
    #=
    prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    for i in 1:Int64(reals.N_real/2)
        prefft_spec[i] = Float64(imag(reals.val[Int64(reals.N_real/2)+i])/pi)
        prefft_spec[Int64(reals.N_real/2)+i] = Float64(imag(reals.val[i])/pi)
    end

    fft_spec = bfft(prefft_spec)*((reals.omega_max-reals.omega_min)/reals.N_real)
    second_der::Complex{Float64} = 0.0 + 0.0*im

    for i in 2:Int64(reals.N_real/2)
        t = 2*pi*(i-1)/(reals.omega_max-reals.omega_min)
        second_der += t^4 * fft_spec[i] * fft_spec[reals.N_real-i+2] /(reals.omega_max-reals.omega_min)
    end
    
    second_der = second_der*2

    lambda::Float64 = 1e-6
    func::Complex{Float64} = abs(1-fft_spec[1])^2 + lambda*second_der
    =#
    
    lambda::Float64 = 1e-6
    tot_int::Float64 = 0.0
    second_square::Float64 = 0.0
    for i in 1:reals.N_real
        tot_int += imag(reals.val[i]/pi)*((reals.omega_max-reals.omega_min)/reals.N_real)
        second_der_i = (imag(reals.val[mod1(i+1,reals.N_real)]/pi) + imag(reals.val[mod1(i-1,reals.N_real)]/pi) - 2*imag(reals.val[i]/pi))*((reals.omega_max-reals.omega_min)/reals.N_real)^(-2)
        second_square += (second_der_i)^2*((reals.omega_max-reals.omega_min)/reals.N_real)
    end
    
    func::Float64 = (1-tot_int) + lambda*second_square
    
    #=
    prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    for i in 1:Int64(reals.N_real/2)
        prefft_spec[i] = Float64(imag(reals.val[Int64(reals.N_real/2)+i])/pi)
        prefft_spec[Int64(reals.N_real/2)+i] = Float64(imag(reals.val[i])/pi)
    end

    fft_spec = bfft(prefft_spec)*2.0*reals.omega_max/reals.N_real

    second_der_prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    for i in 1:Int64(reals.N_real/2)
        t = pi*(i-1)/reals.omega_max
        second_der_prefft_spec[i] = fft_spec[i]*t^2
        second_der_prefft_spec[reals.N_real-i+1] = fft_spec[reals.N_real-i+1]*t^2
    end

    second_der_fft_spec = fft(second_der_prefft_spec)/(2.0*reals.omega_max)

    second_der::Complex{Float64} = 0.0 + 0.0*im

    for i in 1:reals.N_real
        #second_der += second_der_fft_spec[i]*second_der_fft_spec[i]*2.0*reals.omega_max/reals.N_real
        second_der += abs(second_der_fft_spec[i])*2.0*reals.omega_max/reals.N_real
    end
    
    lambda::Float64 = 1e-4
    func::Float64 = abs(1-fft_spec[1])^2 + lambda*second_der
    =#
end


function calc_functional_naive(reals::RealDomainData, abcd::Array{Complex{BigFloat},3}, H::Int64, ab_coeff::Vector{Float64} )

    function calc_g(i)
        z::Complex{BigFloat} = reals.freq[i]
        param::Complex{BigFloat} = 0.0+0.0*im
        for k in 1:H
            param += ab_coeff[k]*Nevanlinna.hardy_basis(z,k-1)
            param += ab_coeff[k+H]*conj(Nevanlinna.hardy_basis(z,k-1))
        end

        theta = (abcd[1,1,i]*param + abcd[1,2,i]) / (abcd[2,1,i]*param + abcd[2,2,i])
        green = im * (one(BigFloat) + theta) / (one(BigFloat) - theta)
    end

    A = collect(imag(calc_g(i))/pi for i in 1:reals.N_real)

    second_square = 0.0
    tot_int = sum(A)*((2.0*reals.omega_max)/reals.N_real)

    #prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    #for i in 1:Int64(reals.N_real/2)
    #    prefft_spec[i] = Float64(A[Int64(reals.N_real/2)+i])
    #    prefft_spec[Int64(reals.N_real/2)+i] = Float64(A[i])
    #end

    #fft_spec = bfft(prefft_spec)*2.0*reals.omega_max/reals.N_real
    fft_spec = bfft(Float64.(A))*2.0*reals.omega_max/reals.N_real

    function calc_second(i)
        (2*pi*(i-1)/(2.0*reals.omega_max))^4 * abs(fft_spec[i])^2 /(2.0*reals.omega_max)
    end

    second_der = sum(collect(calc_second(i) for i in 2:Int64(reals.N_real/2)))

    second_der = second_der*2

    lambda::Float64 = 1e-5
    func::Complex{Float64} = abs(1-tot_int)^2 + lambda*second_der

end


function evaluation(imag::ImagDomainData, reals::RealDomainData, phis::Vector{Complex{BigFloat}}, H::Int64, ab_coeff::Vector{Complex{BigFloat}})
    for i in 1:reals.N_real
        result = Matrix{Complex{BigFloat}}(I, 2, 2) 
        z::Complex{BigFloat} = reals.freq[i]
        for j in 1:imag.N_imag
            prod = Array{Complex{BigFloat}}(undef, 2, 2)
            prod[1,1] = (z - imag.freq[j]) / (z - conj(imag.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imag.freq[j]) / (z - conj(imag.freq[j]))
            prod[2,2] = one(BigFloat)
            result *= prod
        end
        
        param::Complex{BigFloat} = 0.0+0.0*im
        for k in 1:H
            param += ab_coeff[k]*hardy_basis(z,k-1)
            param += ab_coeff[k+H]*conj(hardy_basis(z,k-1))
        end
        
        theta = (result[1,1]*param + result[1,2]) / (result[2,1]*param + result[2,2])
        reals.val[i] = im * (one(BigFloat) + theta) / (one(BigFloat) - theta)
    end
end

function calc_functional(reals::RealDomainData, abcd::Array{Complex{BigFloat},3}, H::Int64, ab_coeff::Vector{Complex{BigFloat}}, hardy_matrix::Array{Complex{BigFloat},2})::Float64
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    green = im * (one(BigFloat) .+ theta) ./ (one(BigFloat) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = sum(A)*((2.0*reals.omega_max)/reals.N_real)

    #本来はこうすべきだが、以下のようにしても結果は変わらない
    #prefft_spec = Array{Complex{Float64}}(undef, reals.N_real)
    #for i in 1:Int64(reals.N_real/2)
    #    prefft_spec[i] = Float64(A[Int64(reals.N_real/2)+i])
    #    prefft_spec[Int64(reals.N_real/2)+i] = Float64(A[i])
    #end
    #fft_spec = bfft(prefft_spec)*2.0*reals.omega_max/reals.N_real
 
    fft_spec = bfft(A)*2.0*reals.omega_max/reals.N_real

    preder_spec = fft_spec[1:Int64(reals.N_real/2)]
    t_vec = 2*pi*Vector(0:Int64(reals.N_real/2)-1)/(2.0*reals.omega_max)

    second_der = 2*sum(t_vec.^4 .* abs.(preder_spec).^2 /(2*reals.omega_max))

    lambda::Float64 = 1e-5
    func::Float64 = abs(1-tot_int)^2 + lambda*second_der

    return func
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

    """
    tot_int = sum(A)*((2.0*reals.omega_max)/reals.N_real)

    fft_spec = bfft(A)*2.0*reals.omega_max/reals.N_real

    preder_spec = fft_spec[1:Int64(reals.N_real/2)]
    t_vec = 2*pi*Vector(0:Int64(reals.N_real/2)-1)/(2.0*reals.omega_max)

    second_der = 2*sum(t_vec.^4 .* abs.(preder_spec).^2 /(2*reals.omega_max))
    """
    tot_int = integrate(reals.freq, A)
    second_der = integrate_squared_second_deriv(reals.freq, A) 

    max_theta = findmax(abs.(param))[1]
    func = abs(1-tot_int)^2 + lambda*second_der

    return func
end


function Nevanlinna_Schur(N_imag::Int64, 
                    omega::Array{T,1}, 
                    green::Array{Complex{T},1},
                    N_real::Int64,
                    omega_max::Float64,
                    eta::Float64,
                    ab_coeff::Array{ComplexF64,1},
                    H::Int64,
                    iter_tol::Int64,
                    lambda::Float64,
                    verbose::Bool=false)::Tuple{ImagDomainData{T}, RealDomainData{T}, Array{ComplexF64,1}, Bool, Bool} where {T<:Real}
    if N_real%2 == 1
        error("N_real must be even number!")
    end
    
    imags = ImagDomainData(N_imag, omega, green)
    reals = RealDomainData(N_real, omega_max, eta, T=T)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    hardy_matrix = calc_hardy_matrix(reals, H)
    
    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(reals, abcd, H, x, hardy_matrix, lambda=lambda)
    end
    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    #function functional(x::Vector{Complex{T}})::Float64
    #    return calc_functional(reals, abcd, H, x, hardy_matrix)
    #end
    #function jacobian(J::Vector{Complex{T}}, x::Vector{Complex{T}})
    #    J .= gradient(functional, x)[1] 
    #end
   
    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = verbose))
    
    if  !(Optim.converged(res))
        #error("Faild to optimize!")
        println("Faild to optimize!")
    end
    
    causality = evaluation(reals, abcd, H, Optim.minimizer(res), hardy_matrix)
    
    return imags, reals, Optim.minimizer(res), causality, (Optim.converged(res))
end

"""
function Nevanlinna_Schur(N_imag::Int64, 
                    omega::Array{T,1}, 
                    green::Array{Complex{T},1},
                    N_real::Int64,
                    omega_max::Float64,
                    eta::Float64,
                    H::Int64;
                    verbose::Bool=false,
                    iterations = 100000,
                    lambda = 1e-5
                    )::Tuple{ImagDomainData{T}, RealDomainData{T}} where {T<:Real}
    if N_real%2 == 1
        error("N_real must be even number!")
    end
    
    imags = ImagDomainData(N_imag, omega, green)
    reals = RealDomainData(N_real, omega_max, eta, T=T)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)
    hardy_matrix = calc_hardy_matrix(reals, H)
    
    ab_coeff  = zeros(ComplexF64, 2*H) 
    
    function functional(x::Vector{ComplexF64})::Float64
        return Nevanlinna.calc_functional(reals, abcd, H, x, hardy_matrix, lambda=lambda)
    end
    
    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end
   
    if verbose
        res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                        Optim.Options(iterations = iterations,
                                      show_trace = true))
    else 
        res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                        Optim.Options(iterations = iterations,
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
"""


