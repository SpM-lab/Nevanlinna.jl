function calc_functional(reals::Nevanlinna.RealDomainData{T},
                         abcd::Array{Complex{T},3},
                         H::Int64,
                         ab_coeff::Vector{S},
                         hardy_matrix::Array{Complex{T},2},
                         lambda::Float64
                         )::Float64 where {S<:Real, T<:Real}
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    green = im * (one(T) .+ theta) ./ (one(T) .- theta)
    A = Float64.(imag(green)./pi)

    second_der = Nevanlinna.integrate_squared_second_deriv(reals.freq, A)

    func = lambda*second_der
    
    return func
end

function functional_opt(ab_coeff::Vector{S},
                        grad    ::Vector{S},
                        reals::Nevanlinna.RealDomainData{T},
                        abcd::Array{Complex{T},3},
                        H::Int64,
                        hardy_matrix::Array{Complex{T},2},
                        lambda::Float64
                        )::Float64 where {S<:Real, T<:Real}
    if length(grad) > 0
        grad .= gradient(x->calc_functional(reals, abcd, H, x, hardy_matrix, lambda), ab_coeff)[1]
    end
    
    return calc_functional(reals, abcd, H, ab_coeff, hardy_matrix, lambda)
    
end
function calc_tot(reals::RealDomainData{T},
                  abcd::Array{Complex{T},3},
                  H::Int64,
                  ab_coeff::Vector{S},
                  hardy_matrix::Array{Complex{T},2}
                  )::Float64 where {S<:Real, T<:Real}
    param = hardy_matrix*ab_coeff

    theta = (abcd[1,1,:].* param .+ abcd[1,2,:]) ./ (abcd[2,1,:].*param .+ abcd[2,2,:])
    green = im * (one(T) .+ theta) ./ (one(T) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = Nevanlinna.integrate(reals.freq, A)

    return (tot_int - 1.0)
end

function tot_opt(ab_coeff::Vector{S},
                 grad    ::Vector{S},
                 reals::Nevanlinna.RealDomainData{T},
                 abcd::Array{Complex{T},3},
                 H::Int64,
                 hardy_matrix::Array{Complex{T},2}
                 )::Float64 where {S<:Real, T<:Real}
    if length(grad) > 0
        grad .= gradient(x->calc_tot(reals, abcd, H, x, hardy_matrix), ab_coeff)[1]
    end

    return calc_tot(reals, abcd, H, ab_coeff, hardy_matrix)
end

function calc_max_theta(reals::RealDomainData{T},
                        abcd::Array{Complex{T},3},
                        H::Int64,
                        ab_coeff::Vector{S},
                        hardy_matrix::Array{Complex{T},2},
                        )::Float64 where {S<:Real, T<:Real}
    param = hardy_matrix*ab_coeff

    max_theta = Float64(findmax(abs.(param))[1])

    return (max_theta - (1.0-1e-10))
end

function max_theta_opt(ab_coeff::Vector{S},
                       grad    ::Vector{S},
                       reals::Nevanlinna.RealDomainData{T},
                       abcd::Array{Complex{T},3},
                       H::Int64,
                       hardy_matrix::Array{Complex{T},2}
                       )::Float64 where {S<:Real, T<:Real}
    if length(grad) > 0
        grad .= gradient(x->calc_max_theta(reals, abcd, H, x, hardy_matrix), ab_coeff)[1]
    end

    return calc_max_theta(reals, abcd, H, ab_coeff, hardy_matrix)
end
