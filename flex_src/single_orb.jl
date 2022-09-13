function int_n(
        mu   ::Float64, 
        sekf ::Array{ComplexF64,3}, 
        lat  ::SingleLatticeModel,
        basis::FiniteTempBasisSet
        )::Float64 

    new_gkf = make_int_giw(sekf, mu, lat, basis)

    g_f = dropdims(sum(new_gkf, dims=(2,3)), dims=(2,3)) / (lat.nsize^2)

    g_l = fit(basis.smpl_wn_f, g_f, dim=1)
    g_tau0 = dot(basis.basis_f.u(0), g_l)

    n  = 1.0 + real(g_tau0)
    n  = 2.0 * n #for spin

    return real(n)
end

"""
Calculate the chemical potential for interacting system
"""
function calc_int_chem(
        sekf ::Array{ComplexF64,3}, 
        lat  ::SingleLatticeModel, 
        basis::FiniteTempBasisSet;
        minmu = -5.0,
        maxmu = 5.0
        )::Float64

    f = x -> int_n(x, sekf, lat, basis) - lat.filling
    return Roots.find_zero(f, (minmu, maxmu), Roots.Brent())
end

"""
Make free green function on fermionic sampling frequencies
"""
function make_free_giw(
        lat  ::SingleLatticeModel, 
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,3}

    return make_free_giw(lat, SparseIR.β(basis), basis.smpl_wn_f.sampling_points)
end

"""
Make free green function on given fermionic frequencies
"""
function make_free_giw(
        lat    ::SingleLatticeModel, 
        beta   ::Float64,
        vsample::Vector{FermionicFreq}
        )::Array{ComplexF64,3}

    giw = Array{ComplexF64}(undef, length(vsample), lat.nsize, lat.nsize)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:length(vsample)
        iv = valueim(vsample[iw], beta)
        locg = iv - lat.hk[ix, iy] + lat.freemu
        giw[iw,ix,iy] = 1.0/locg 
    end
    giw
end

"""
Make interacting green function for IR basis
"""
function make_int_giw(
        sekf ::Array{ComplexF64,3}, 
        mu   ::Float64,
        lat  ::SingleLatticeModel,
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,3}

    giw = Array{ComplexF64}(undef, getfnw(basis), lat.nsize, lat.nsize)
    beta = SparseIR.β(basis)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:getfnw(basis)
        iv::ComplexF64 = valueim(basis.smpl_wn_f.sampling_points[iw], beta)
        locg = iv - lat.hk[ix, iy] + mu - sekf[iw,ix,iy]
        giw[iw, ix, iy] = 1.0/locg
    end
    giw
end


"""
find and return appropriate U where alpha_s < 1. (initU < U < maxU)
"""
function find_next_U(
        current_solution::SingleSCFSolution,
        try_ratio_U     ::Float64,
        lat             ::SingleLatticeModel,
        basis           ::FiniteTempBasisSet,
        max_depth=10
        )::Float64

    if max_depth > 0 && calc_stoner_with_green(current_solution.gkf, try_ratio_U, lat, basis) > (1.0 - current_solution.precision)
        #Try again with narrow range since maxU is too large
        finalU = try_ratio_U*lat.U
        curU   = current_solution.ratio_U*lat.U
        try_ratio_U = 0.5*(finalU+curU)/lat.U
        return find_next_U(current_solution, try_ratio_U, lat, basis, max_depth-1)
    end
    try_ratio_U
end

function flex_scheme(
        gkf      ::Array{ComplexF64,3}, 
        ratio_U  ::Float64,
        precision::Float64,
        lat      ::SingleLatticeModel,
        var      ::Variables,
        basis    ::FiniteTempBasisSet 
        )::Tuple{Variables, Array{ComplexF64,3}, Float64}

    # Compute chi0 
    chi0kf = calc_chi0(gkf, lat, basis)

    # Compute stoner factor
    stoner = calc_stoner_with_chi0(chi0kf, ratio_U, lat, basis)
    println("alpha_s=$(stoner)")

    chiskf = chi0kf ./ (1.0 .- (ratio_U*lat.U) .* chi0kf) #spin susceptibility
    chickf = chi0kf ./ (1.0 .+ (ratio_U*lat.U) .* chi0kf) #charge susceptibility

    #Compute effective interaction for particle-hole
    Vkf = calc_Vph(chiskf, chickf, chi0kf, ratio_U, lat, basis)

    # Compute self-energy using sparse modeling
    sekf = calc_se(gkf, Vkf, lat, basis)

    #determination of shift of chemical potential
    mu_new = calc_int_chem(sekf, lat, basis)
    println("chemical potential shift=$(mu_new-lat.freemu)")
    new_gkf_tmp = make_int_giw(sekf, mu_new, lat, basis)
    
    #impose G^\daggar(iomega,k) = G(-iomega,k)
    new_gkf = Array{ComplexF64}(undef, getfnw(basis), lat.nsize, lat.nsize)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:getfnw(basis)
        new_gkf[iw,ix,iy] = 0.5 * (new_gkf_tmp[iw,ix,iy] + conj(new_gkf_tmp[getfnw(basis)-iw+1,ix,iy])) 
    end

    #evaluate error between previous and new green function
    max = findmax(abs.(gkf.-new_gkf))[1]
    println("max=$(max), U=$(ratio_U*lat.U)")

    return_gkf = var.mixing .* new_gkf .+ (1.0-var.mixing) .* gkf

    new_var = Variables(var.mixing, max, stoner, var.prestoner)

    return new_var, return_gkf, mu_new
end

"""
Compute chi0 
"""
function calc_chi0(
        gkf  ::Array{ComplexF64,3}, 
        lat  ::SingleLatticeModel,
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,3}

    @assert size(gkf) == (getfnw(basis), lat.nsize, lat.nsize)
    ntau_ = getntau(basis)

    # G(iw, k) -> G(l, k)
    gkl = fit(basis.smpl_wn_f, gkf, dim=1)
    @assert size(gkl) == (getnl(basis), lat.nsize, lat.nsize)

    # G(l, k) -> G(l, r)
    grl = ifft(gkl, [2,3])

    # G(l, r) -> G(tau, r)
    grt::Array{ComplexF64,3} = evaluate(basis.smpl_tau_f, grl, dim=1)
    @assert size(grt) == (getntau(basis), lat.nsize, lat.nsize)

    # Make irreducible susceptibility
    chi0rt = Array{ComplexF64}(undef, ntau_, lat.nsize, lat.nsize)

    # Since nsize-ix+2 over nsize when ix or iy == 1,
    # we enter mod1 function. mod1(nsize+1,nsize)=1 and mod(ix ,nsize)=ix (for 2<=ix<=nsize) !
    for iy in 1:lat.nsize, ix in 1:lat.nsize, itau in 1:ntau_
        chi0rt[itau, ix, iy] = grt[itau, ix, iy] * grt[ntau_-itau+1, mod1(lat.nsize-ix+2,lat.nsize), mod1(lat.nsize-iy+2,lat.nsize)]
    end

    #chi0(tau,r) -> chi0(l,r)
    chi0rl = fit(basis.smpl_tau_b, chi0rt, dim=1)
    @assert size(chi0rl) == (getnl(basis), lat.nsize, lat.nsize)

    #chi0(l,r) -> chi0(l,k)
    chi0kl = fft(chi0rl, [2,3])

    #chi0(l,k) -> chi0(iv,k)
    chi0kf = evaluate(basis.smpl_wn_b, chi0kl, dim=1)
    @assert size(chi0kf) == (getbnw(basis), lat.nsize, lat.nsize)

    return chi0kf
end


"""
Compute effective interaction for particle-hole channel
"""
function calc_Vph(
        chiskf ::Array{ComplexF64,3}, 
        chickf ::Array{ComplexF64,3}, 
        chi0kf ::Array{ComplexF64,3}, 
        ratio_U::Float64,
        lat    ::SingleLatticeModel,
        basis  ::FiniteTempBasisSet
        )::Array{ComplexF64,3}

    V = Array{ComplexF64}(undef, getbnw(basis), lat.nsize, lat.nsize)
    V .= (ratio_U*lat.U)^2 .* (1.5.*chiskf .+ 0.5.*chickf .- chi0kf)

    return V
end

"""
Compute self-energy 
"""
function calc_se(
        gkf  ::Array{ComplexF64,3},
        Vkf  ::Array{ComplexF64,3}, 
        lat  ::SingleLatticeModel,
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,3}

    @assert size(gkf) == (getfnw(basis), lat.nsize, lat.nsize)

    # G(iw, k) -> G(l, k)
    gkl = fit(basis.smpl_wn_f, gkf, dim=1)
    @assert size(gkl) == (getnl(basis), lat.nsize, lat.nsize)

    # G(l, k) -> G(l, r)
    grl = ifft(gkl, [2,3])

    # G(l, r) -> G(tau, r)
    grt = evaluate(basis.smpl_tau_f, grl, dim=1)
    @assert size(grt) == (getntau(basis), lat.nsize, lat.nsize)
    @assert size(Vkf) == (getbnw(basis), lat.nsize, lat.nsize)

    # V(iv, k) -> V(l, k)
    Vkl = fit(basis.smpl_wn_b, Vkf, dim=1)
    @assert size(Vkl) == (getnl(basis), lat.nsize, lat.nsize)

    # V(l, k) -> V(l, r)
    Vrl = ifft(Vkl, [2,3])

    # Warning!! 
    #Now, tau sampling points for fermion and boson is same. There is the possibility of the change of implementation. Therefore, we must evaluate the sampling points with fermionic sampling points explicitly!! 
    # V(l, r) -> V(tau, r)
    loc_smpl_tau_b = TauSampling(basis.basis_b, basis.smpl_tau_f.sampling_points)
    Vrt = evaluate(loc_smpl_tau_b, Vrl, dim=1)
    @assert size(Vrt) == (length(loc_smpl_tau_b.sampling_points), lat.nsize, lat.nsize)

    #make self-energy
    sert = zeros(ComplexF64, getntau(basis), lat.nsize, lat.nsize)

    for iy in 1:lat.nsize, ix in 1:lat.nsize, itau in 1:getntau(basis)
        sert[itau, ix, iy] += Vrt[itau, ix, iy] * grt[itau, ix, iy]
    end

    #se(tau,r) -> se(l,r)
    serl = fit(basis.smpl_tau_f, sert, dim=1)
    @assert size(serl) == (getnl(basis), lat.nsize, lat.nsize)

    #se(l,r) -> se(l,k)
    sekl = fft(serl, [2,3]) 

    #se(l,k) -> se(iw,k)
    sekf = evaluate(basis.smpl_wn_f, sekl, dim=1)
    @assert size(sekf) == (getfnw(basis), lat.nsize, lat.nsize)

    return sekf
end

"""
Calc stoner factor
"""
function calc_stoner_with_green(
        gkf    ::Array{ComplexF64,3},
        ratio_U::Float64,
        lat    ::SingleLatticeModel,
        basis  ::FiniteTempBasisSet
        )::Float64

    # Compute Chi0 
    chi0kf = calc_chi0(gkf, lat, basis)

    return calc_stoner_with_chi0(chi0kf, ratio_U, lat, basis)
end

function calc_stoner_with_chi0(
        chi0   ::Array{ComplexF64,3},
        ratio_U::Float64,
        lat    ::SingleLatticeModel,
        basis  ::FiniteTempBasisSet
        )::Float64

    chis = chi0 ./ (1.0 .- ratio_U*lat.U*chi0) #spin susceptibility
    stoner = findmax(abs.(ratio_U*lat.U*chi0))[1]

    return stoner
end

"""
Compute self-consistent solutions for next_U
return converged, solution
"""
function scf_next_U(
        sol         ::SingleSCFSolution,
        next_ratio_U::Float64,
        lat         ::SingleLatticeModel,
        var         ::Variables,
        basis       ::FiniteTempBasisSet,
        max_loop    ::Int64,
        )::Tuple{Variables, Bool, SingleSCFSolution}

    if sol.pre_calc
        loc_precision = sol.loose_precision
    else
        loc_precision = sol.precision
    end

    gkf_copy = copy(sol.gkf)
    mu_new::Float64 = lat.freemu
    converged::Bool = false

    for iloop in 1:max_loop
        var_new, gkf_new, mu_new = flex_scheme(gkf_copy, next_ratio_U, loc_precision, lat, var, basis)
        gkf_copy = gkf_new
        var = var_new
        if var.stoner > 1
           error("Stoner factor become greater than unity. Please reduce mixing parameter!")
        end
        if var.max < loc_precision
            # Converged!
            converged = true
            break
        end
    end

    if !converged
        error("Not converged!")
    end

    return var, converged, SingleSCFSolution(sol.pre_calc, sol.precision, sol.loose_precision, next_ratio_U, mu_new, gkf_copy)
end

function flex_exe(
        sol  ::SingleSCFSolution,
        lat  ::SingleLatticeModel,
        var  ::Variables,
        basis::FiniteTempBasisSet
        )::Tuple{Variables,SingleSCFSolution}
    
    while sol.ratio_U*lat.U < lat.U
        println("Current value of U = $(sol.ratio_U*lat.U)")
        
        #determine next_U
        next_ratio_U = find_next_U(sol, 1.0, lat, basis)
        println("Trying U = $(next_ratio_U*lat.U) with mixing = $(var.mixing)")

        #Loosely Compute self-consistent solution for next_U 
        var_new, converged, sol_new = scf_next_U(sol, next_ratio_U, lat, var, basis, 3000)
        var = var_new
    
        #update mixing parameter
        var.mixing = new_mixing(var.mixing, var.prestoner, var.stoner)
        var.prestoner = var.stoner
    
        sol = sol_new
    end
    
    #Fully compute self-consistent solutions for final_U
    full_sol = SingleSCFSolution(false, sol.precision, sol.loose_precision, 1.0, sol.mu, sol.gkf)
    println("full calc U = $(full_sol.ratio_U*lat.U) with mixing = $(var.mixing)")
    var_new, converged, full_sol = scf_next_U(full_sol, 1.0, lat, var, basis, 3000)
    var = var_new
    
    if !converged
        error("Not converged!")
    end
    
    return var, full_sol
end


"""
Compute chi0 (only for tests)
"""
function calc_chi0_brute_force(beta::Float64, latt::SingleLatticeModel, nmatsu::Int64)
    nsize = latt.nsize
    nallsize = nsize^2
    T = 1/beta
    vsample = FermionicFreq.(2*collect(0:nmatsu-1) .+ 1)

    #make green function
    giw = make_free_giw(latt, beta, vsample)

    chi0 = zeros(Complex{Float64}, (nsize, nsize))
    for iy1=1:nsize, ix1=1:nsize
        for iy2=1:nsize, ix2=1:nsize, iw=1:nmatsu
            #正のomegaと負のomegaに対する和で2倍 G*(iw,k) = -G(-iw,k)
            chi0[ix1,iy1] += -2.0*real(giw[iw,mod1(-ix1+ix2+1,nsize),mod1(-iy1+iy2+1,nsize)] * giw[iw,ix2,iy2]) * (T/nallsize)
        end
    end
    return chi0
end
