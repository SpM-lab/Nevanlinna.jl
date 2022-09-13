abstract type AbstractDecomposed end
struct Decomposed1P <: AbstractDecomposed
    freqidof::Array{ComplexF64,4}
    mom::Array{ComplexF64,3}
    function Decomposed1P(freqidof, mom)
        size(freqidof)[end] == size(mom)[end] || error("Bond dimension mismatch")
        new(freqidof, mom)
    end
    function Decomposed1P(ntau, nidof, nsize, D)
        new(
            Array{ComplexF64}(undef, ntau, nidof, nidof, D),
            Array{ComplexF64}(undef, nsize, nsize, D)
        )
    end
end

struct Decomposed2P <: AbstractDecomposed
    freqidof::Array{ComplexF64,6}
    mom::Array{ComplexF64,3}
    function Decomposed2P(freqidof, mom)
        size(freqidof)[end] == size(mom)[end] || error("Bond dimension mismatch")
        new(freqidof, mom)
    end
    function Decomposed2P(ntau, nidof, nsize, D)
        new(
            Array{ComplexF64}(undef, ntau, nidof, nidof, nidof, nidof, D),
            Array{ComplexF64}(undef, nsize, nsize, D)
        )
    end
end

getn(obj::AbstractDecomposed) = size(obj.freqidof, 1)
getnidof(obj::AbstractDecomposed) = size(obj.freqidof, 2)
getnsize(obj::AbstractDecomposed) = size(obj.mom, 1)
getbonddim(obj::AbstractDecomposed) = size(obj.freqidof)[end]

Base.copy(x::Decomposed2P) = Decomposed2P(copy(x.freqidof), copy(x.mom))
Base.copy(x::Decomposed1P) = Decomposed1P(copy(x.freqidof), copy(x.mom))

function compress(x::Decomposed2P; eps=1e-14)
    D = getbonddim(x)
    n, nidof, nsize = getn(x), getnidof(x), getnsize(x)

    mps = MPS{ComplexF64}(
        [reshape(x.freqidof, :, 1, D), reshape(x.mom, :, D, 1)]
        )
    compress!(mps, eps)
    distributecoeff!(mps)
    new_D = rightbonddim(mps, 1)
    return Decomposed2P(
        reshape(mps.tensors[1], n, nidof, nidof, nidof, nidof, new_D),
        reshape(mps.tensors[2], nsize, nsize, new_D),
    )
end


function _full(x::AbstractDecomposed)
    D = getbonddim(x)
    size1 = size(x.freqidof)[1:end-1]
    size2 = size(x.mom)[1:end-1]
    x1 = reshape(x.freqidof, (:, D))
    x2 = reshape(x.mom, (:, D))
    return reshape(x1 * transpose(x2), size1..., size2...)
end

function _to_mps(x::T) where {T<:AbstractDecomposed}
    D = getbonddim(x)
    return MPS{ComplexF64}([reshape(x.freqidof, :, 1, D), reshape(x.mom, :, D, 1)])
end

norm2(x::AbstractDecomposed) = norm2(_to_mps(x))

full(x::Decomposed1P) =
    permutedims(_full(x), (1, 4, 5, 2, 3))

full(x::Decomposed2P) =
    permutedims(_full(x), (1, 6, 7, 2, 3, 4, 5))


function _calc_chi0rt_decomp!(chi0rt::Decomposed2P, grt_decomp::Decomposed1P)
    ntau = getn(grt_decomp)
    nidof = getnidof(grt_decomp)
    nsize = getnsize(grt_decomp)
    D = getbonddim(grt_decomp)
    revr = mod1.((nsize-ix+2 for ix in 1:nsize), nsize)

    grt_decomp_right1 = grt_decomp.freqidof[end:-1:1, :, :, :]
    grt_decomp_right2 = grt_decomp.mom[revr, revr, :]
    grt_decomp1, grt_decomp2 = grt_decomp.freqidof, grt_decomp.mom

    chi0rt1 = reshape(
        chi0rt.freqidof, ntau, nidof, nidof, nidof, nidof, D, D)
    chi0rt2 = reshape(
        chi0rt.mom, nsize, nsize, D, D)

    for id in 1:nidof, ic in 1:nidof, ib in 1:nidof, ia in 1:nidof, τ in 1:ntau
        for e in 1:D, d in 1:D
            chi0rt1[τ, ia, ib, ic, id, d, e] = grt_decomp1[τ, ia, ic, d] * grt_decomp_right1[τ, id, ib, e]
        end
    end

    for y in 1:nsize, x in 1:nsize
        for e in 1:D, d in 1:D
            chi0rt2[x, y, d, e] = grt_decomp2[x, y, d] * grt_decomp_right2[x, y, e]
        end
    end
    nothing
end


function calc_chi0_decomp(
        gkf::Array{ComplexF64,5}, 
        lat::LatticeModel,
        basis::FiniteTempBasisSet;
        eps=1e-16
        )::Decomposed2P
    ntau = getntau(basis)
    nidof = lat.nidof
    nsize = lat.nsize
    fnw = getfnw(basis)

    @assert size(gkf) == (fnw, nsize, nsize, nidof, nidof)

    # gkf_decomp: (fnw, nidof, nidof, D), (nsize, nsize, D)
    gkf_factors = decomp(
        reshape(
            permutedims(gkf, (1, 4, 5, 2, 3)), (getfnw(basis) * nidof^2, nsize^2)),
        eps
    )
    gkf_decomp = Decomposed1P(
        reshape(gkf_factors[1], fnw, nidof, nidof, :),
        reshape(gkf_factors[2], (nsize, nsize, :))
    )
    D = getbonddim(gkf_decomp)

    # G(iw, k) -> G(l, k)
    gkl = Decomposed1P(
            fit(basis.smpl_wn_f, gkf_decomp.freqidof, dim=1),
            gkf_decomp.mom)

    # G(l, k) -> G(l, r) -> G(tau, r)
    grt_decomp = Decomposed1P(
            evaluate(basis.smpl_tau_f, gkl.freqidof, dim=1),
            ifft(gkl.mom, [1,2]))

    # make irreducible susceptibility
    # chi0rt: (ntau, nidof, nidof, nidof, nidof, D^2), (nsize, nsize, D^2)
    chi0rt = Decomposed2P(ntau, nidof, nsize, D^2)
    _calc_chi0rt_decomp!(chi0rt, grt_decomp)

    # Compression
    chi0rt = compress(chi0rt)
    #println("debugl ", getbonddim(chi0rt), " ", getbonddim(chi0rt_new))

    # chi0(tau,r) -> chi0(l,r) -> chi0(l,k)
    # chi0kl: (nl, nidof, nidof, nidof, nidof, D^2), (nsize, nsize, D^2)
    chi0kl = Decomposed2P(
        fit(basis.smpl_tau_b, chi0rt.freqidof, dim=1),
        fft(chi0rt.mom, [1,2])
    )

    # chi0(l,k) -> chi0(iv,k)
    # chi0kf: (fnw, nidof, nidof, nidof, nidof, D^2), (nsize, nsize, D^2)
    chi0kf = Decomposed2P(
        evaluate(basis.smpl_wn_b, chi0kl.freqidof, dim=1),
        chi0kl.mom
    )

    return chi0kf
end


function _mul(a::Decomposed2P, b::Decomposed2P)
    Da = getbonddim(a)
    Db = getbonddim(b)
    bnw = getn(a)
    nsize = getnsize(a)
    nidof = getnidof(a)
    nidof2 = nidof^2

    mom_ab = reshape(ein"xyd,xye->xyde"(a.mom, b.mom), nsize, nsize, Da * Db)

    # (-1, -2, 1, -4) * (-1, 1, -3, -5) = (-1, -2, -3,  -4, -5)
    # (w, a, b, D) * (w, b, c, E) = (w, a, c, D, E)
    freqidof_ab = zeros(ComplexF64, bnw, nidof2, nidof2, Da, Db)
    freqidof_a = reshape(a.freqidof, (bnw, nidof2, nidof2, Da))
    freqidof_b = reshape(b.freqidof, (bnw, nidof2, nidof2, Db))
    for iw in 1:bnw, da in 1:Da, db in 1:Db
        freqidof_ab[iw, :, :, da, db] =
            view(freqidof_a, iw, :, :, da) *
            view(freqidof_b, iw, :, :, db)
    end
    freqidof_ab = reshape(freqidof_ab, (bnw, nidof, nidof, nidof, nidof, Da * Db))

    return Decomposed2P(freqidof_ab, mom_ab)
end

function Base.:+(a::Decomposed2P, b::Decomposed2P)
    # (w, a, b, c, d, D)
    freqidof = cat(a.freqidof, b.freqidof, dims=6)
    # (kx, ky, D)
    mom = cat(a.mom, b.mom, dims=3)
    return Decomposed2P(freqidof, mom)
end

function Base.:-(a::Decomposed2P, b::Decomposed2P)
    # (w, a, b, c, d, D)
    freqidof = cat(a.freqidof, -b.freqidof, dims=6)
    # (kx, ky, D)
    mom = cat(a.mom, b.mom, dims=3)
    return Decomposed2P(freqidof, mom)
end


function Base.:*(a::Number, b::Decomposed2P)
    return Decomposed2P(sqrt(a) * b.freqidof, sqrt(a) * b.mom)
end

"""
Compute χ0 + χ0 U χ
"""
function _update_chi(
        chi0kf::Decomposed2P,
        chikf::Decomposed2P,
        ratio_U::Float64,
        lat::LatticeModel,
        basis::FiniteTempBasisSet)::Decomposed2P
    bnw = getbnw(basis)
    nidof = lat.nidof
    D = getbonddim(chikf)
    nidof2 = nidof^2

    # Compute U * chi
    freqidof = reshape(chikf.freqidof, bnw, nidof2, nidof2, D)
    freqidof = permutedims(freqidof, (2, 3, 1, 4)) # (nidof2, nidof2, bnw, D)
    freqidof = reshape(freqidof, (nidof2, nidof2 * bnw * D))
    freqidof_Uχ = reshape((ratio_U * lat.matU) * freqidof, (nidof, nidof, nidof, nidof, bnw, D))
    freqidof_Uχ = permutedims(freqidof_Uχ, (5, 1, 2, 3, 4, 6))
    Uχ = Decomposed2P(freqidof_Uχ, chikf.mom)

    return chi0kf + _mul(chi0kf, Uχ)
end

"""
Compute chi iteratively
"""
function calc_chi(
        chi0kf::Decomposed2P,
        chikf_ini::Decomposed2P,
        ratio_U::Float64,
        lat::LatticeModel,
        basis::FiniteTempBasisSet;
        eps_compress = 1e-10,
        tol = 1e-10, # tol for terminating iterations
        maxiter = 1000
        )::Decomposed2P
    chikf = copy(chikf_ini)
    chikf_prev = chikf_ini
    for ite in 1:maxiter
        chikf = compress(_update_chi(chi0kf, chikf, ratio_U, lat, basis), eps=eps_compress)
        if norm2(chikf - chikf_prev) < tol * norm2(chikf)
            break
        end
        chikf_prev = chikf
    end
    return chikf
end