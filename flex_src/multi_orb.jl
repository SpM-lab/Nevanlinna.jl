"""
Compute total density of electrons for given mu, assuming PM
"""
#=
function free_n(
        mu::Float64,
        mesh::Int64,
        lat::LatticeModel,
        basis::FiniteTempBasisSet)::Float64

    N::Float64 = 0
    H = Array{ComplexF64,2}(undef, lat.nidof, lat.nidof)
    for ix in 1:mesh, iy in 1:mesh
        kx::Float64 = (2*π*(ix-1))/mesh
        ky::Float64 = (2*π*(iy-1))/mesh
        H .= Hami(kx, ky, lat.nidof, lat.t_pra, lat.alpha)
        for i in 1:lat.nidof
            H[i,i] -= mu
        end
        diagH = eigen!(H)
        N += sum(fermi_dirac.(real.(diagH.values), basis.beta))
    end
    n = N/(mesh*mesh)
    return n
end
=#

function int_n(
        mu   ::Float64, 
        sekf ::Array{ComplexF64,5}, 
        lat  ::LatticeModel,
        basis::FiniteTempBasisSet
        )::Float64 

    new_gkf = make_int_giw(sekf, mu, lat, basis)
    g_f = zeros(ComplexF64, size(new_gkf, 1))
    for ikx in 1:lat.nsize, iky in 1:lat.nsize, i in 1:lat.nidof
        g_f .+= @view new_gkf[:, ikx, iky, i, i]
    end
    g_tau0 = dot(
        basis.basis_f.u(0),
        fit(basis.smpl_wn_f, g_f, dim=1)
    ) / lat.nsize^2
    return lat.nidof + real(g_tau0)
end

"""
Calculate the chemical potential for interacting system
"""
function calc_int_chem(
        sekf ::Array{ComplexF64,5}, 
        lat  ::LatticeModel, 
        basis::FiniteTempBasisSet;
        minmu = -5.0,
        maxmu = 5.0
        )::Float64

    f  = x -> int_n(x, sekf, lat, basis) - lat.filling
    return Roots.find_zero(f, (minmu, maxmu), Roots.Brent())
end

"""
Make free green function for IR basis
"""
function make_free_giw(
        lat  ::LatticeModel, 
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,5}

    return make_free_giw(lat, SparseIR.β(basis), basis.smpl_wn_f.sampling_points)
end

"""
Make free green function on given frequencies
"""
function make_free_giw(
        lat    ::LatticeModel, 
        beta   ::Float64,
        vsample::Vector{FermionicFreq}
        )::Array{ComplexF64,5}

    giw = Array{ComplexF64}(undef, length(vsample), lat.nsize, lat.nsize, lat.nidof, lat.nidof)
    locg = Array{ComplexF64}(undef, lat.nidof, lat.nidof)
    work = InvWorkspace(locg)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:length(vsample)
        iv::ComplexF64 = valueim(vsample[iw], beta)
        for j in 1:lat.nidof, i in 1:lat.nidof
            locg[i, j] = -lat.hk[ix, iy, i, j]
        end
        for i in 1:lat.nidof
            locg[i, i] += iv + lat.freemu
        end
        inv!(locg, work)
        for j in 1:lat.nidof, i in 1:lat.nidof
            giw[iw, ix, iy, i, j] = locg[i, j]
        end
    end
    return giw
end

"""
Make interacting green function for IR basis
"""
function make_int_giw(
        sekf ::Array{ComplexF64,5}, 
        mu   ::Float64,
        lat  ::LatticeModel,
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,5}

    giw = Array{ComplexF64}(undef, getfnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)
    locg = Array{ComplexF64}(undef, lat.nidof, lat.nidof)
    work = InvWorkspace(locg)
    beta = SparseIR.β(basis)
    fnw = getfnw(basis)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:fnw
        iv::ComplexF64 = valueim(basis.smpl_wn_f.sampling_points[iw], beta)
        for j in 1:lat.nidof, i in 1:lat.nidof
            locg[i, j] = -lat.hk[ix, iy, i, j] - sekf[iw, ix, iy, i, j]
        end
        for i in 1:lat.nidof
            locg[i, i] += iv + mu 
        end
        inv!(locg, work)
        for j in 1:lat.nidof, i in 1:lat.nidof
            giw[iw, ix, iy, i, j] = locg[i, j]
        end
    end
    giw
end


"""
find and return appropriate U where alpha_s < 1. (initU < U < maxU)
"""
function find_next_U(
        current_solution::SCFSolution,
        try_ratio_U     ::Float64,
        lat             ::LatticeModel,
        basis           ::FiniteTempBasisSet,
        max_depth=10)::Float64

    if max_depth > 0 && calc_stoner_with_green(current_solution.gkf, try_ratio_U, lat, basis) > (1.0 - current_solution.precision)
        #Try again with narrow range since maxU is too large
        finalU = try_ratio_U
        curU   = current_solution.ratio_U
        try_ratio_U = 0.5*(finalU+curU)
        return find_next_U(current_solution, try_ratio_U, lat, basis, max_depth-1)
    end
    try_ratio_U
end

function flex_scheme(
        gkf      ::Array{ComplexF64,5}, 
        ratio_U  ::Float64,
        precision::Float64,
        lat      ::LatticeModel,
        var      ::Variables,
        basis    ::FiniteTempBasisSet 
        )::Tuple{Variables, Array{ComplexF64,5}, Float64}

    # Compute chi0 
    chi0kf = calc_chi0(gkf, lat, basis)

    # Compute stoner factor
    stoner = calc_stoner_with_chi0(chi0kf, ratio_U, lat, basis)
    println("alpha_s=$(stoner)")

    # Compute chi
    chikf = calc_chi(chi0kf, ratio_U, lat, basis)

    #Compute effective interaction for particle-hole
    Vkf = calc_Vph(chikf, chi0kf, ratio_U, lat, basis)

    # Compute self-energy using sparse modeling
    sekf = calc_se(gkf, Vkf, lat, basis)

    #determination of shift of chemical potential
    mu_new = calc_int_chem(sekf, lat, basis)
    println("chemical potential shift=$(mu_new-lat.freemu)")
    new_gkf_tmp = make_int_giw(sekf, mu_new, lat, basis)
    
    #impose G^\daggar(iomega,k) = G(-iomega,k)
    new_gkf = Array{ComplexF64}(undef, getfnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:getfnw(basis)
        new_gkf[iw,ix,iy,:,:] .= (view(new_gkf_tmp,iw,ix,iy,:,:) .+ (view(new_gkf_tmp,getfnw(basis)-iw+1,ix,iy,:,:))') .* 0.5
    end

    #evaluate error between previous and new green function
    max = findmax(abs.(gkf.-new_gkf))[1]
    println("max=$(max), ratio_U=$(ratio_U)")

    return_gkf = var.mixing .* new_gkf .+ (1.0-var.mixing) .* gkf

    new_var = Variables(var.mixing, max, stoner, var.prestoner)

    return new_var, return_gkf, mu_new
end

"""
Compute chi0 
"""
function calc_chi0(
        gkf  ::Array{ComplexF64,5}, 
        lat  ::LatticeModel,
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,7}

    ntau = getntau(basis)

    @assert size(gkf) == (getfnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    # G(iw, k) -> G(l, k)
    gkl = fit(basis.smpl_wn_f, gkf, dim=1)
    @assert size(gkl) == (getnl(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    # G(l, k) -> G(l, r)
    grl = ifft(gkl, [2,3])

    # G(l, r) -> G(tau, r)
    grt::Array{ComplexF64,5} = evaluate(basis.smpl_tau_f, grl, dim=1)
    @assert size(grt) == (ntau, lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    #make irreducible susceptibility
    chi0rt = Array{ComplexF64}(undef, ntau, lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)
    _calc_chi0rt!(chi0rt, grt)

    #chi0(tau,r) -> chi0(l,r)
    chi0rl = fit(basis.smpl_tau_b, chi0rt, dim=1)
    @assert size(chi0rl) == (getnl(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)

    #chi0(l,r) -> chi0(l,k)
    chi0kl = fft(chi0rl, [2,3])

    #chi0(l,k) -> chi0(iv,k)
    chi0kf = evaluate(basis.smpl_wn_b, chi0kl, dim=1)
    @assert size(chi0kf) == (getbnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)

    return chi0kf
end

"""
Compute chi0 
"""
function _calc_chi0rt!(
    chi0rt::Array{ComplexF64,7},
    grt::Array{ComplexF64,5})
    ntau = size(grt, 1)
    nsize = size(grt, 2)
    nidof = size(grt, 4)
    for id in 1:nidof, ic in 1:nidof, ib in 1:nidof, ia in 1:nidof
        for iy in 1:nsize, ix in 1:nsize, itau in 1:ntau
            chi0rt[itau, ix, iy, ia, ib, ic, id] = grt[itau, ix, iy, ia, ic] * grt[ntau-itau+1, mod1(nsize-ix+2, nsize), mod1(nsize-iy+2, nsize), id, ib]
        end
    end
end

"""
Compute χ0 + χ0 U χ
"""
function _update_chi(
        chi0kf::Array{ComplexF64,7},
        chikf::Array{ComplexF64,7},
        ratio_U::Float64,
        lat    ::LatticeModel,
        basis  ::FiniteTempBasisSet
        )::Array{ComplexF64,7}

    nsize = lat.nsize
    nidof = lat.nidof
    nidof2 = nidof^2

    res = Array{ComplexF64}(undef, getbnw(basis), nsize, nsize, nidof, nidof, nidof, nidof)
    loc_chi0 = Array{ComplexF64}(undef, nidof2, nidof2)
    loc_chi  = Array{ComplexF64}(undef, nidof2, nidof2)
    loc_res  = Array{ComplexF64}(undef, nidof2, nidof2)
    loc_res  = Array{ComplexF64}(undef, nidof2, nidof2)
    for iy in 1:nsize, ix in 1:nsize, iw in 1:getbnw(basis)
        loc_chi0 .= reshape(
            view(chi0kf,iw,ix,iy,:,:,:,:), nidof2, nidof2)
        loc_chi .= reshape(
            view(chikf,iw,ix,iy,:,:,:,:), nidof2, nidof2)
        loc_res .= loc_chi0 * (ratio_U*lat.matU) * loc_chi
        res[iw,ix,iy,:,:,:,:] .= reshape(loc_res, nidof, nidof, nidof, nidof)
    end
    return chi0kf + res
end


"""
Compute chi
"""
function calc_chi(
        chi0kf::Array{ComplexF64,7},
        ratio_U::Float64,
        lat::LatticeModel,
        basis::FiniteTempBasisSet)::Array{ComplexF64,7}
    nsize = lat.nsize
    nidof = lat.nidof
    nidof = lat.nidof
    nidof2 = nidof^2

    chikf = Array{ComplexF64}(undef, getbnw(basis), nsize, nsize, nidof, nidof, nidof, nidof)
    loc_chi0 = Array{ComplexF64}(undef, nidof2, nidof2)
    loc_chi  = Array{ComplexF64}(undef, nidof2, nidof2)
    inv_stn = Array{ComplexF64}(undef, nidof2, nidof2)
    matunit = Matrix{ComplexF64}(I, nidof2, nidof2)
    work = InvWorkspace(matunit)
    for iy in 1:nsize, ix in 1:nsize, iw in 1:getbnw(basis)
        loc_chi0 .= reshape(view(chi0kf,iw,ix,iy,:,:,:,:), nidof2, nidof2)
        mul!(inv_stn, loc_chi0, ratio_U*lat.matU)
        inv_stn .= matunit .- inv_stn
        inv!(inv_stn, work)
        mul!(loc_chi, inv_stn, loc_chi0)
        chikf[iw,ix,iy,:,:,:,:] .= reshape(loc_chi, nidof, nidof, nidof, nidof)
    end
    return chikf
end


"""
Compute effective interaction for particle-hole channel
"""
function calc_Vph(
        chikf  ::Array{ComplexF64,7}, 
        chi0kf ::Array{ComplexF64,7}, 
        ratio_U::Float64,
        lat    ::LatticeModel,
        basis  ::FiniteTempBasisSet
        )::Array{ComplexF64,7}

    tmp1 = Array{ComplexF64}(undef, lat.nidof*lat.nidof, lat.nidof*lat.nidof)
    tmp2 = Array{ComplexF64}(undef, lat.nidof*lat.nidof, lat.nidof*lat.nidof)
    V = Array{ComplexF64}(undef, getbnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:getbnw(basis)
        tmp1 .= reshape(view(chikf,iw,ix,iy,:,:,:,:) , lat.nidof*lat.nidof, lat.nidof*lat.nidof) 
        tmp1 .-= 0.5.* reshape(view(chi0kf,iw,ix,iy,:,:,:,:), lat.nidof*lat.nidof, lat.nidof*lat.nidof) 
        mul!(tmp2, tmp1, ratio_U*lat.matU)
        mul!(tmp1, ratio_U*lat.matU, tmp2)
        V[iw,ix,iy,:,:,:,:] .= reshape(tmp1, lat.nidof, lat.nidof, lat.nidof, lat.nidof)
    end
    
    return V
end

"""
Compute self-energy 
"""
function calc_se(
        gkf  ::Array{ComplexF64,5},
        Vkf  ::Array{ComplexF64,7}, 
        lat  ::LatticeModel,
        basis::FiniteTempBasisSet
        )::Array{ComplexF64,5}

    ntau = getntau(basis)
    @assert size(gkf) == (getfnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    # G(iw, k) -> G(l, k)
    gkl = fit(basis.smpl_wn_f, gkf, dim=1)
    @assert size(gkl) == (getnl(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    # G(l, k) -> G(l, r)
    grl = ifft(gkl, [2,3])

    # G(l, r) -> G(tau, r)
    grt = evaluate(basis.smpl_tau_f, grl, dim=1)
    @assert size(grt) == (ntau, lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    @assert size(Vkf) == (getbnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)

    # V(iv, k) -> V(l, k)
    Vkl = fit(basis.smpl_wn_b, Vkf, dim=1)
    @assert size(Vkl) == (getnl(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)

    # V(l, k) -> V(l, r)
    Vrl = ifft(Vkl, [2,3])

    # V(l, r) -> V(tau, r)
    Vrt = evaluate(basis.smpl_tau_b, Vrl, dim=1)
    @assert size(Vrt) == (ntau, lat.nsize, lat.nsize, lat.nidof, lat.nidof, lat.nidof, lat.nidof)

    #make self-energy
    sert = zeros(ComplexF64, ntau, lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    for ib in 1:lat.nidof, ia in 1:lat.nidof
        for id in 1:lat.nidof, ic in 1:lat.nidof
          for iy in 1:lat.nsize, ix in 1:lat.nsize, itau in 1:ntau
            sert[itau, ix, iy, ia, ib] += Vrt[itau, ix, iy, ia, ic, ib, id] * grt[itau, ix, iy, ic, id]
          end
        end
    end

    #se(tau,r) -> se(l,r)
    serl = fit(basis.smpl_tau_f, sert, dim=1)
    @assert size(serl) == (getnl(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    #se(l,r) -> se(l,k)
    sekl = fft(serl, [2,3]) 

    #se(l,k) -> se(iw,k)
    sekf = evaluate(basis.smpl_wn_f, sekl, dim=1)
    @assert size(sekf) == (getfnw(basis), lat.nsize, lat.nsize, lat.nidof, lat.nidof)

    return sekf
end

"""
Calc stoner factor
"""
function calc_stoner_with_green(
        gkf    ::Array{ComplexF64,5},
        ratio_U::Float64,
        lat    ::LatticeModel,
        basis  ::FiniteTempBasisSet
        )::Float64

    # Compute Chi0 
    chi0kf = calc_chi0(gkf, lat, basis)

    return calc_stoner_with_chi0(chi0kf, ratio_U, lat, basis)
end

function calc_stoner_with_chi0(
        chi0   ::Array{ComplexF64,7},
        ratio_U::Float64,
        lat    ::LatticeModel,
        basis  ::FiniteTempBasisSet
        )::Float64

    loc_chi0 = Array{ComplexF64,2}(undef, lat.nidof*lat.nidof, lat.nidof*lat.nidof)
    tmp = Array{ComplexF64,2}(undef, lat.nidof*lat.nidof, lat.nidof*lat.nidof)
    val = Vector{ComplexF64}(undef,lat.nidof*lat.nidof)
    lwork = lwork_check(tmp, val)
    W = EigWorkspace(tmp, lwork)
    mat_stoner = Array{ComplexF64}(undef, getbnw(basis), lat.nsize, lat.nsize, lat.nidof*lat.nidof)
    for iy in 1:lat.nsize, ix in 1:lat.nsize, iw in 1:getbnw(basis)
        loc_chi0 .= reshape(view(chi0,iw,ix,iy,:,:,:,:), lat.nidof*lat.nidof, lat.nidof*lat.nidof)
        mul!(tmp, loc_chi0, ratio_U*lat.matU)
        eig!(tmp, val, W)
        mat_stoner[iw,ix,iy,:] .= val 
    end
    stoner = findmax(abs.(mat_stoner))[1]

    return stoner
end

"""
Compute self-consistent solutions for next_U
return converged, solution
"""
function scf_next_U(
        sol         ::SCFSolution,
        next_ratio_U::Float64,
        lat         ::LatticeModel,
        var         ::Variables,
        basis       ::FiniteTempBasisSet,
        max_loop    ::Int64,
        )::Tuple{Variables, Bool, SCFSolution}

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
    return var, converged, SCFSolution(sol.pre_calc, sol.precision, sol.loose_precision, next_ratio_U, mu_new, gkf_copy)
end

function flex_exe(
        sol  ::SCFSolution,
        lat  ::LatticeModel,
        var  ::Variables,
        basis::FiniteTempBasisSet
        )::Tuple{Variables,SCFSolution}
    
    while sol.ratio_U < 1.0
        println("Current value of raito_U = $(sol.ratio_U)")
        
        #determine next_U
        next_ratio_U = find_next_U(sol, 1.0, lat, basis)
        println("Trying ratio_U = $(next_ratio_U) with mixing = $(var.mixing)")

        #Loosely Compute self-consistent solution for next_U 
        var_new, converged, sol_new = scf_next_U(sol, next_ratio_U, lat, var, basis, 3000)
        var = var_new
    
        #update mixing parameter
        var.mixing = new_mixing(var.mixing, var.prestoner, var.stoner)
        var.prestoner = var.stoner
    
        sol = sol_new
    end
    
    #Fully compute self-consistent solutions for final_U
    full_sol = SCFSolution(false, sol.precision, sol.loose_precision, 1.0, sol.mu, sol.gkf)
    println("full calc raito_U = $(full_sol.ratio_U) with mixing = $(var.mixing)")
    var_new, converged, full_sol = scf_next_U(full_sol, 1.0, lat, var, basis, 3000)
    var = var_new
    
    if !converged
        error("Not converged!")
    end
    
    return var, full_sol
end

