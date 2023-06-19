module Nevanlinna                                                                                                   
using LinearAlgebra
using GenericLinearAlgebra
using Optim
using Zygote
using Comonicon
using TOML

# Some hack
using MultiFloats
function Base.convert(t::Type{Float64x2}, ::Irrational{:π})
    return Float64x2(BigFloat(π, precision=128))
end

function Float64x2(::Irrational{:π})
    return Float64x2(BigFloat(π, precision=128))
end

#==
# log2 is not implemented in MultiFloats.
# This will use the BigFloat implementation of
# log2, which will not be as fast as a pure-MultiFloat implementation.
MultiFloats.use_bigfloat_transcendentals()
==#

include("export.jl")
include("util.jl")
include("data.jl")
include("solver.jl")
include("ham_solver.jl")
include("hardy.jl")
#include("nevanlinna_impl.jl")
include("core.jl")
include("schur.jl")
include("optimize.jl")


@cast function bare(input_data::String, input_param::String, output_data::String)
    f= open(input_data, "r")
    data = split.(readlines(f),'\t')
    close(f)
    wn = im.*parse.(BigFloat, collect(Iterators.flatten(data))[1:3:end])
    gw = parse.(BigFloat,collect(Iterators.flatten(data))[2:3:end]) .+ im.*parse.(BigFloat,collect(Iterators.flatten(data))[3:3:end])

    param = TOML.parsefile(input_param)

    N_real::Int64     = param["basic"]["N_real"]
    w_max::Float64    = param["basic"]["w_max"]
    eta::Float64      = param["basic"]["eta"]
    sum_rule::Float64 = param["basic"]["sum_rule"]
    H_max::Int64      = param["basic"]["H_max"]
    iter_tol::Int64   = param["basic"]["iter_tol"]
    lambda::Float64   = param["basic"]["lambda"]

    verbose::Bool      = false
    pick_check::Bool   = true
    optimization::Bool = true
    ini_iter_tol::Int64 = 500
    mesh::Symbol       = :linear

    if haskey(param, "option")
        if haskey(param["option"], "verbose")
            verbose = param["option"]["verbose"]
        end
        if haskey(param["option"], "pick_check")
            pick_check = param["option"]["pick_check"]
        end
        if haskey(param["option"], "optimization")
            optimization = param["option"]["optimization"]
        end
        if haskey(param["option"], "ini_iter_tol")
            optimization = param["option"]["ini_iter_tol"]
        end
        if haskey(param["option"], "mesh")
            mesh = param["option"]["mesh"]
        end
    end

    sol = NevanlinnaSolver(wn, gw, N_real, w_max, eta, sum_rule, H_max, iter_tol, lambda, verbose=verbose, pick_check = pick_check, optimization=optimization, ini_iter_tol=ini_iter_tol, mesh=mesh)
    if optimization 
        solve!(sol)
    end

    open(output_data, "w") do f
        for i in 1:N_real
            println(f, Float64(real(sol.reals.freq[i])), "\t", Float64(imag(sol.reals.val[i]))/pi)
        end
    end
end

@cast function hamburger(input_data::String, input_moment::String, input_param::String, output_data::String)
    f= open(input_data, "r")
    data = split.(readlines(f),'\t')
    close(f)
    wn = im.*parse.(BigFloat, collect(Iterators.flatten(data))[1:3:end])
    gw = parse.(BigFloat,collect(Iterators.flatten(data))[2:3:end]) .+ im.*parse.(BigFloat,collect(Iterators.flatten(data))[3:3:end])

    f = open(input_moment, "r")
    moments = Complex{BigFloat}.(parse.(Float64, readlines(f)))
    close(f)

    param = TOML.parsefile(input_param)

    N_real::Int64     = param["basic"]["N_real"]
    w_max::Float64    = param["basic"]["w_max"]
    eta::Float64      = param["basic"]["eta"]
    sum_rule::Float64 = param["basic"]["sum_rule"]
    H_max::Int64      = param["basic"]["H_max"]
    iter_tol::Int64   = param["basic"]["iter_tol"]
    lambda::Float64   = param["basic"]["lambda"]

    verbose::Bool      = false
    pick_check::Bool   = true
    optimization::Bool = true
    mesh::Symbol       = :linear

    if haskey(param, "option")
        if haskey(param["option"], "verbose")
            verbose = param["option"]["verbose"]
        end
        if haskey(param["option"], "pick_check")
            pick_check = param["option"]["pick_check"]
        end
        if haskey(param["option"], "optimization")
            optimization = param["option"]["optimization"]
        end
        if haskey(param["option"], "ini_iter_tol")
            optimization = param["option"]["ini_iter_tol"]
        end
        if haskey(param["option"], "mesh")
            mesh = param["option"]["mesh"]
        end
    end

    sol = HamburgerNevanlinnaSolver(moments, wn, gw, N_real, w_max, eta, sum_rule, H_max, iter_tol, lambda, verbose=verbose, pick_check = pick_check, optimization = optimization, ini_iter_tol=ini_iter_tol, mesh = mesh)
    if optimization 
        solve!(sol)
    end

    open(output_data, "w") do f
        for i in 1:N_real
            println(f, Float64(real(sol.nev_st.reals.freq[i])), "\t", Float64(imag(sol.val[i]))/pi)
        end
    end
end

@main
end
