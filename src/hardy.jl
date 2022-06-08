function hardy_basis(z::Complex{BigFloat}, k::Int64)
    #(z-im)^k/(sqrt(pi)*(z+im)^(k+1))
    w = (z-im)/(z+im)
    0.5*im*(w^(k+1)-w^k)/(sqrt(pi))
end

function hardy_basis(x::Float64, y::Float64, k::Int64)
    z::Complex{BigFloat} = big(x) +im*big(y)
    
    return hardy_basis(z, k)
end
