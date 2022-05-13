using SparseIR
using PyPlot
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
rcParams["font.size"] = 16
rcParams["text.latex.preamble"] = raw"\usepackage{amsmath}"
using LinearAlgebra
using Revise
using Optim
using FFTW
using ForwardDiff
include("./Nevanlinna.jl")

# Three Gaussian peaks (normalized to 1)
gaussian(x, mu, sigma) = exp(-((x-mu)/sigma)^2)/(sqrt(π)*sigma)

#rho(omega) = 0.2*gaussian(omega, 0.0, 0.15) +
#    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

rho(omega) = gaussian(omega, 0.0, 0.15)
#    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

beta = 100
wmax = 1000
IR_basis_set = FiniteTempBasisSet(beta, wmax, 1e-7)
#basis = FiniteTempBasis(fermion, beta, wmax, 1e-7)

rhol = [overlap(IR_basis_set.basis_f.v[l], rho) for l in 1:size(IR_basis_set.basis_f)]
gl = - IR_basis_set.basis_f.s .* rhol

gw = evaluate(IR_basis_set.smpl_wn_f, gl)
hnw = Int64(length(IR_basis_set.smpl_wn_f.sampling_points)/2)

setprecision(128)

test_gw   = Array{Complex{BigFloat}}(undef, hnw)
test_smpl = Array{BigFloat}(undef, hnw)

for i in 1:hnw
    test_smpl[i]= parse(BigFloat, string(IR_basis_set.smpl_wn_f.sampling_points[hnw+i]*pi/beta))
    test_gw[i]  = parse(BigFloat, string(real(gw[hnw+i]))) + parse(BigFloat, string(imag(gw[hnw+i])))*im
end

N_imag    =  12
N_real    =  6000
omega_min = -10.0
omega_max =  10.0
eta       =  0.001
H         =  25
ab_coeff = zeros(Float64, 2*H)
imaginary = ImagDomainData(N_imag, test_smpl, test_gw)
phis = calc_phis(imaginary)
reals = RealDomainData(N_real, omega_min, omega_max, eta)
evaluation(imaginary, reals, phis, H, ab_coeff)

function functional(imag::ImagDomainData, reals::RealDomainData, phis::Vector{Complex{BigFloat}}, H::Int64, ab_coeff::Vector{Float64})
    evaluation(imaginary, reals, phis, H, ab_coeff)
    return real(calc_functional(reals))
end


ff = x->functional(imaginary, reals, phis, H, x)

@time res = optimize(ff, ab_coeff, BFGS())

println(res)

param_f = open( "optim_param.dat", "w")
for i in 1:2*H
    println(param_f, "$(Optim.minimizer(res)[i])")
end
close(param_f)
