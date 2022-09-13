module IRFLEX

using FFTW
import LinearAlgebra: svd, eigen!, I, mul!, dot, qr
using SparseArrays
using SparseIR
using Roots
using TensorOperations
import SparseIR: value, valueim
import OMEinsum: @ein, @ein_str

function __init__()
    # Some initialization
end

#include("mps.jl")
include("struct.jl")
include("single_struct.jl")
include("exports.jl")
include("basis.jl")
include("Hamiltonian.jl")
include("util.jl")
include("single_orb.jl")
include("multi_orb.jl")
include("multi_orb_decomp.jl")
#include("irbasis.jl")
#include("multi_orb_old.jl")
#include("multi_orb_new.jl")

end
