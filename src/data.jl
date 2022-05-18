struct ImagDomainData
    N_imag::Int64
    freq  ::Array{Complex{BigFloat},1}
    val   ::Array{Complex{BigFloat},1}
end

function ImagDomainData(N_imag::Int64,
                        matsu ::Array{BigFloat,1},
                        green ::Array{Complex{BigFloat},1}
                        )::ImagDomainData
    val  = Array{Complex{BigFloat}}(undef, N_imag) 
    freq = Array{Complex{BigFloat}}(undef, N_imag) 
    
    for i in 1:N_imag
        freq[i] = matsu[i]*im
        val[i]  = (-green[i] - im) / (-green[i] + im)
    end
    
    Pick = Array{Complex{BigFloat}}(undef, N_imag, N_imag)
    
    for j in 1:N_imag
        for i in 1:N_imag
            freq_i = (freq[i] - im) / (freq[i] + im)
            freq_j = (freq[j] - im) / (freq[j] + im)
            nom = one(BigFloat) - val[i] * conj(val[j])
            den = one(BigFloat) - freq_i * conj(freq_j)
            Pick[i,j] = nom / den
        end
        Pick[j,j] += big(1e-250)
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

struct RealDomainData
    N_real   ::Int64
    omega_max::Float64
    eta      ::Float64
    freq     ::Array{Complex{BigFloat},1}
    val      ::Array{Complex{BigFloat},1}
end

function RealDomainData(N_real   ::Int64,
                        omega_max::Float64,
                        eta      ::Float64
                        )::RealDomainData
    val  = Array{Complex{BigFloat}}(undef, N_real) 
    freq = Array{Complex{BigFloat}}(undef, N_real) 
    
    inter::BigFloat = parse(BigFloat, string((2.0*omega_max) / (N_real-1)))
    temp ::BigFloat = parse(BigFloat, string(-omega_max))
    
    freq[1] = parse(BigFloat, string(-omega_max)) + parse(BigFloat, string(eta))*im
    for i in 2:N_real
        temp += inter
        freq[i] = temp + parse(BigFloat, string(eta))*im
    end
    
    return RealDomainData(N_real, omega_max, eta, freq, val)
end
