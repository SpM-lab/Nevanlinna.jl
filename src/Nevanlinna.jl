module Nevanlinna                                                                                                   
using LinearAlgebra
#using FFTW
using SparseIR
using Optim
using Zygote

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
include("hardy.jl")
include("nevanlinna_impl.jl")
include("schur.jl")

end
