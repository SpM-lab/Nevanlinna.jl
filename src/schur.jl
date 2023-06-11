function calc_opt_N_imag(N      ::Int64,
                         wn     ::Array{Complex{T},1},
                         gw     ::Array{Complex{T},1};
                         verbose::Bool=false
                         )::Int64 where {T<:Real}
    @assert N == length(wn)
    @assert N == length(gw)

    freq = (wn  .- im) ./ (wn  .+ im)
    val  = (-gw .- im) ./ (-gw .+ im)

    k::Int64 = 0
    success::Bool = true

    while success
        k += 1
        Pick = Array{Complex{T}}(undef, k, k)

        for j in 1:k
            for i in 1:k
                num = one(T) - val[i]  * conj(val[j])
                den = one(T) - freq[i] * conj(freq[j])
                Pick[i,j] = num / den
            end

            Pick[j,j] += T(1e-250)
        end

        success = issuccess(cholesky(Pick,check = false))

        if k == N
            break
        end
    end

    if verbose
        if !(success)
            println("N_imag is setted as $(k-1)")
        else
            println("N_imag is setted as $(N)")
        end
    end

    if !(success)
        return (k-1)
    else
        return (N)
    end
end


