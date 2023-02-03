struct ImagDomainData{T<:Real}
    N_imag::Int64               #The number of points used in Nevanlinna algorithm
    freq  ::Array{Complex{T},1} #The values of Matsubara frequencies
    val   ::Array{Complex{T},1} #The values of negative of Green function
end

function ImagDomainData(wn     ::Array{Complex{T},1},
                        gw     ::Array{Complex{T},1},
                        N_imag ::Int64;
                        verbose::Bool = false
                        )::ImagDomainData{T} where {T<:Real}

    val  = Array{Complex{T}}(undef, N_imag) 
    freq = Array{Complex{T}}(undef, N_imag) 
    
    for i in 1:N_imag
        freq[i] = wn[i]
        val[i]  = (-gw[i] - im) / (-gw[i] + im) 
    end
    
    Pick = Array{Complex{T}}(undef, N_imag, N_imag)
    
    for j in 1:N_imag
        for i in 1:N_imag
            freq_i = (freq[i] - im) / (freq[i] + im)
            freq_j = (freq[j] - im) / (freq[j] + im)
            nom = one(T) - val[i] * conj(val[j])
            den = one(T) - freq_i * conj(freq_j)
            Pick[i,j] = nom / den
        end
        Pick[j,j] += T(1e-250)
    end
    
    success = issuccess(cholesky(Pick,check = false))
    
    if verbose
        if success
            println("Pick matrix is positive semi-definite.")
        else
            println("Pick matrix is non positive semi-definite matrix in Schur method.")
        end
    end
    
    freq = reverse(freq)
    val  = reverse(val)
    
    return ImagDomainData(N_imag, freq, val)
end

struct RealDomainData{T<:Real}
    N_real  ::Int64               #The number of mesh in real axis
    w_max   ::Float64             #The energy cutoff of real axis
    eta     ::Float64             #The paramer. The retarded Green function is evaluated at omega+i*eta
    sum_rule::Float64             #The value of sum of spectral function
    freq    ::Array{Complex{T},1} #The values of frequencies of retarded Green function
    val     ::Array{Complex{T},1} #The values of negative of retarded Green function
end

function RealDomainData(N_real  ::Int64,
                        w_max   ::Float64,
                        eta     ::Float64,
                        sum_rule::Float64
                        ;
                        T::Type=BigFloat,
                        small_omega::Float64 = 1e-5,
                        mesh::Symbol=:linear
                        )::RealDomainData{T}

    if mesh === :linear
        val = Array{Complex{T}}(collect(LinRange(-w_max, w_max, N_real)))
        freq = val .+ eta * im
        return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
    elseif mesh === :log
        half_N = N_real ÷ 2
        mesh = exp.(LinRange(log.(small_omega), log.(w_max), half_N))
        val = Array{Complex{T}}([reverse(-mesh); mesh])
        freq = val .+ eta * im
        return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
    elseif mesh === :test
        val  = Array{Complex{T}}(undef, N_real) 
        freq = Array{Complex{T}}(undef, N_real) 
        inter::T = big(2.0*w_max) / (N_real-1)
        temp ::T = big(-w_max)
        freq[1] = -big(w_max) + big(eta)*im
        for i in 2:N_real
            temp += inter
            freq[i] = temp + big(eta)*im
        end
        return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
    else
        throw(ArgumentError("Invalid mesh"))
    end
end
