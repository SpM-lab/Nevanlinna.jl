struct SingleLatticeModel
    nsize  ::Int64
    filling::Float64
    freemu ::Float64             #chemical potential for free system
    U      ::Float64
    hk     ::Array{ComplexF64,2} #pre-calculated Hamiltonian

    function SingleLatticeModel(
            nsize  ::Int64,
            filling::Float64,
            U      ::Float64,
            beta   ::Float64,
            hami   ::Function
            )::SingleLatticeModel
    
        #Compute Hamiltonian
        hk = Array{ComplexF64,2}(undef, nsize, nsize)
        for iy in 1:nsize, ix in 1:nsize
            kx::Float64 = (2*π*(ix-1))/nsize
            ky::Float64 = (2*π*(iy-1))/nsize
            hk[ix, iy] = hami(kx, ky)
        end
        
        #Compute free chemical potential
        freemu = calc_free_chem(filling, beta, hami)
    
        new(nsize, filling, freemu, U, hk)
    end
end

"""
self-consistent solution at ratio_U*U
"""
struct SingleSCFSolution
    pre_calc       ::Bool                       #pre_calc=true  -> loose    pre_calc=false -> full
    precision      ::Float64                    #precision for full-calculation
    loose_precision::Float64                    #precision for loose-calculation
    ratio_U        ::Float64                    #ratio of interaction between this solution and final_U
    mu             ::Float64                    #chemical potential for interacting system
    gkf            ::Array{Complex{Float64},3}  #Green function
end

function free_n(
        mu  ::Float64, 
        mesh::Int64, 
        beta::Float64,
        hami::Function
        )::Float64

    N::Float64 = 0
    for ix in 1:mesh, iy in 1:mesh
        kx::Float64 = (2*π*(ix-1))/mesh
        ky::Float64 = (2*π*(iy-1))/mesh
        e = hami(kx, ky)
        e -= mu
        N += 2.0*fermi_dirac(e, beta)
    end
    n = N/(mesh*mesh)
    return n
end

"""
Calculate the chemical potential for non-interacting system
"""
function calc_free_chem(
        filling::Float64,
        beta   ::Float64,
        hami   ::Function
        )::Float64

    minmu::Float64 = -5.0
    maxmu::Float64 = 5.0
    mesh_for_chem_search::Int64 = 100
    freemu = find_zero(
        x -> free_n(x, mesh_for_chem_search, beta, hami) - filling,
        minmu, maxmu,
        1e-10,
        10000
    )
    return freemu
end

function single_Hami(x::Float64, y::Float64, t_pra::Float64)::Float64
    e::Float64 = -2.0*(cos(x)+cos(y)) + 4.0*t_pra*cos(x)*cos(y) 
    return e
end
