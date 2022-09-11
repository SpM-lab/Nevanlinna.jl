getfnw(basis::FiniteTempBasisSet) = length(basis.smpl_wn_f.sampling_points)
getbnw(basis::FiniteTempBasisSet) = length(basis.smpl_wn_b.sampling_points)
getntau(basis::FiniteTempBasisSet) = length(basis.smpl_tau_f.sampling_points)
getnl(basis::FiniteTempBasisSet) = length(basis.basis_f)