struct ImagDomainData{T<:Real}
    N_imag::Int64
    freq  ::Array{Complex{T},1}
    val   ::Array{Complex{T},1}
end

function ImagDomainData(N_imag::Int64,
                        matsu ::Array{T,1},
                        green ::Array{Complex{T},1}
                        )::ImagDomainData{T} where {T<:Real}
    val  = Array{Complex{T}}(undef, N_imag) 
    freq = Array{Complex{T}}(undef, N_imag) 
    
    for i in 1:N_imag
        freq[i] = matsu[i]*im
        val[i]  = (-green[i] - im) / (-green[i] + im)
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
    
    if success
        println("Pick matrix is positive semi-definite.")
    else
        println("Pick matrix is non positive semi-definite matrix in Schur method.")
    end
    
    freq = reverse(freq)
    val  = reverse(val)
    
    return ImagDomainData(N_imag, freq, val)
end

struct RealDomainData{T<:Real}
    N_real   ::Int64
    omega_max::Float64
    eta      ::Float64
    freq     ::Array{Complex{T},1}
    val      ::Array{Complex{T},1}
end

function RealDomainData(N_real   ::Int64,
                        omega_max::Float64,
                        eta      ::Float64
                        ;
                        T::Type=BigFloat,
                        small_omega::Float64 = 1e-5,
                        mesh::Symbol=:linear
                        )::RealDomainData{T}
    if mesh === :linear
        val = Array{Complex{T}}(collect(LinRange(-omega_max, omega_max, N_real)))
        freq = val .+ eta * im
        return RealDomainData(N_real, omega_max, eta, freq, val)
    elseif mesh === :log
        half_N = N_real ÷ 2
        mesh = exp.(LinRange(log.(small_omega), log.(omega_max), half_N))
        val = Array{Complex{T}}([reverse(-mesh); mesh])
        freq = val .+ eta * im
        return RealDomainData(N_real, omega_max, eta, freq, val)
    else
        throw(ArgumentError("Invalid mesh"))
    end
    
    """
    val  = Array{Complex{T}}(undef, N_real) 
    freq = Array{Complex{T}}(undef, N_real) 
    inter::T = (2.0*omega_max) / (N_real-1)
    temp ::T = -omega_max
    freq[1] = -omega_max + eta*im
    for i in 2:N_real
        temp += inter
        freq[i] = temp + eta*im
    end
    """
end
