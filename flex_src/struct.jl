struct LatticeModel
    nsize  ::Int64
    nidof  ::Int64
    filling::Float64
    freemu ::Float64                         #chemical potential for free system
    matU   ::SparseMatrixCSC{Float64,Int64}  #matrix for bare interaction
    hk     ::Array{ComplexF64,4}             #pre-calculated Hamiltonian

    # Construct single-band model
    function LatticeModel(
            nsize  ::Int64,
            filling::Float64,
            U      ::Float64,
            beta   ::Float64,
            hami   ::Function
            )::LatticeModel
        nidof = 2

        # Compute bare interaction in sparse matrix
        spI  = [1, 2, 3, 4] #i component of Utensor
        spJ  = [4, 2, 3, 1] #j component of Utensor
        valU = [-U, U, U, -U] # non-zero value of Utensor 
        matU = sparse(spI, spJ, valU)
    
        # Precompute Hamiltonian
        hk = zeros(ComplexF64, nsize, nsize, nidof, nidof)
        for iy in 1:nsize, ix in 1:nsize
            kx::Float64 = (2*π*(ix-1))/nsize
            ky::Float64 = (2*π*(iy-1))/nsize
            hk[ix, iy, :, :] = hami(kx, ky)
        end
        
        # Compute free chemical potential
        freemu = calc_free_chem(nidof, filling, beta, hami)
    
        new(nsize, nidof, filling, freemu, matU, hk)
    end
    
    # Construct general model
    function LatticeModel(
            nsize  ::Int64,
            nidof  ::Int64,
            filling::Float64,
            beta   ::Float64,
            matU   ::SparseMatrixCSC{Float64,Int64},
            hami   ::Function
            )::LatticeModel

        # Precompute Hamiltonian
        hk = zeros(ComplexF64, nsize, nsize, nidof, nidof)
        for iy in 1:nsize, ix in 1:nsize
            kx::Float64 = (2*π*(ix-1))/nsize
            ky::Float64 = (2*π*(iy-1))/nsize
            hk[ix, iy, :, :] = hami(kx, ky)
        end
        
        # Compute free chemical potential
        freemu = calc_free_chem(nidof, filling, beta, hami)
    
        new(nsize, nidof, filling, freemu, matU, hk)
    end

end

mutable struct Variables
    mixing   ::Float64  #ratio of mixing of Green function between new and old one
    max      ::Float64  #maximum value of difference between new and old Green function
    stoner   ::Float64  #stoner factor
    prestoner::Float64  #stoner factor for previous step
end

"""
self-consistent solution at ratio_U*U
"""
struct SCFSolution
    pre_calc       ::Bool                      #pre_calc=true  -> loose    pre_calc=false -> full
    precision      ::Float64                   #precision for full-calculation
    loose_precision::Float64                   #precision for loose-calculation
    ratio_U        ::Float64                   #ratio of interaction between this solution and final_U
    mu             ::Float64                   #chemical potential for interacting system
    gkf            ::Array{Complex{Float64},5} #Green function
end


function free_n(
        mu   ::Float64, 
        mesh ::Int64, 
        nidof::Int64,
        beta ::Float64,
        hami ::Function
        )::Float64

    N::Float64 = 0
    H = Array{ComplexF64,2}(undef, nidof, nidof)
    for ix in 1:mesh, iy in 1:mesh
        kx::Float64 = (2*π*(ix-1))/mesh
        ky::Float64 = (2*π*(iy-1))/mesh
        H .= hami(kx, ky)
        for i in 1:nidof
            H[i,i] -= mu
        end
        diagH = eigen!(H)
        N += sum(fermi_dirac.(real.(diagH.values), beta))
    end
    n = N/(mesh*mesh)
    return n
end

"""
Calculate the chemical potential for non-interacting system
"""
function calc_free_chem(
        nidof  ::Int64,
        filling::Float64,
        beta   ::Float64,
        hami   ::Function;
        minmu  ::Float64 = -5.0,
        maxmu  ::Float64 = 5.0,
        mesh_for_chem_search::Int64 = 100
        )::Float64

    freemu = find_zero(
        x -> free_n(x, mesh_for_chem_search, nidof, beta, hami) - filling,
        minmu, maxmu,
        1e-10,
        10000
    )
    return freemu
end
