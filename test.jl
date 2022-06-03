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
#using ForwardDiff
using Zygote
#using Flux
include("./src/Nevanlinna.jl")

# Three Gaussian peaks (normalized to 1)
gaussian(x, mu, sigma) = exp(-((x-mu)/sigma)^2)/(sqrt(π)*sigma)

#rho(omega) = 0.2*gaussian(omega, 0.0, 0.15) +
#    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

rho(omega) = gaussian(omega, 0.0, 0.15)
#    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

beta = 100
wmax = 1000
IR_basis_set = FiniteTempBasisSet(Float64(beta), Float64(wmax), 1e-7)
#basis = FiniteTempBasis(fermion, beta, wmax, 1e-7)

rhol = [overlap(IR_basis_set.basis_f.v[l], rho) for l in 1:length(IR_basis_set.basis_f)]
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
omega_max =  10.0
eta       =  0.001
H         =  10
ab_coeff = zeros(Float64, 2*H)

imaginary = Nevanlinna.ImagDomainData(N_imag, test_smpl, test_gw)
pre_reals = Nevanlinna.RealDomainData(N_real, omega_max, eta)
reals     = Nevanlinna.RealDomainData(N_real, omega_max, eta)

phis = Nevanlinna.calc_phis(imaginary)
abcd = Nevanlinna.calc_abcd(imaginary, reals, phis)
hardy_matrix = Nevanlinna.calc_hardy_matrix(reals, H)

Nevanlinna.evaluation(imaginary, pre_reals, abcd, H, ab_coeff, hardy_matrix)

f = x->Nevanlinna.calc_functional(reals, abcd, H, x, hardy_matrix)

function j(J::Vector, x)
   J .= gradient(f, x)[1]
end

@time res = optimize(f, j, ab_coeff, BFGS(), Optim.Options(iterations = 100000,
                                                          show_trace = true))

#@time res = optimize(f, j, ab_coeff, ConjugateGradient(), Optim.Options(iterations = 100000,
#                                                          show_trace = true))

println(res)

println(Optim.minimizer(res))







