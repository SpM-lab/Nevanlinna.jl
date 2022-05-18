function hardy_basis(z::Complex{BigFloat}, k::Int64)
    (z-im)^k/(sqrt(pi)*(z+im)^(k+1))
end

function hardy_basis(x::Float64, y::Float64, k::Int64)
    bigx = parse(BigFloat, string(x))
    bigy = parse(BigFloat, string(y))
    
    z::Complex{BigFloat} = x +im*y
    
    return hardy_basis(z, k)
end
