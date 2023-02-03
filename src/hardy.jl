function hardy_basis(z::Complex{T}, k::Int64) where {T<:Real}
    w = (z-im)/(z+im)
    0.5*im*(w^(k+1)-w^k)/(sqrt(pi))
end

#function hardy_basis(x::Float64, y::Float64, k::Int64)
#    z::Complex{T} = T(x) +im*T(y)
#    return hardy_basis(z, k)
#end

function calc_hardy_matrix(reals::RealDomainData{T}, 
                           H::Int64
                           )::Array{Complex{T}, 2} where {T<:Real}
    hardy_matrix = Array{Complex{T}}(undef, reals.N_real, 2*H)
    for k in 1:H
        hardy_matrix[:,2*k-1] .=      hardy_basis.(reals.freq,k-1)
        hardy_matrix[:,2*k]   .= conj(hardy_basis.(reals.freq,k-1))
    end
    return hardy_matrix
end
