function fermi_dirac(x::Float64, beta::Float64)::Float64
    0.5 * (1.0 - tanh(x*beta*0.5))
end

"""
Find a root by the bisection algorithm
The function f(x) has opposite signs at x=xmin and x=xmax, i.e., f(xmin) * f(xmax) < 0.
"""
function find_zero(f::Function, xmin::Float64, xmax::Float64, atol, max_div::Int64=1000)
    success = false
    f_xmin = f(xmin)
    f_xmax = f(xmax)
    if f_xmin * f_xmax >= 0
        error("f(xmin) * f(xmax) >=0 !")
    end
    for i in 1:max_div
        if xmax - xmin < atol
            success = true
            break
        end
        xmid = 0.5 * (xmin + xmax)
        f_xmid = f(xmid)
        if f_xmin * f_xmid > 0
            f_xmin = f_xmid
            xmin = xmid
        else
            f_xmax = f_xmid
            xmax = xmid
        end
    end
    if !success
        error("Faild to find a root!")
    end
    0.5 * (xmin + xmax)
end

"""
Heuristics for new mixing parameter
"""
function new_mixing(mixing::Float64, prestoner::Float64, stoner::Float64)
    r = mixing
    if prestoner < 0.98 && stoner > 0.98
        r = 0.5*r
    elseif prestoner < 0.99 && stoner > 0.99
        r = 0.5*r
    elseif prestoner < 0.995 && stoner > 0.995
        r = 0.5*r
    elseif prestoner < 0.997 && stoner > 0.997
        r = 0.5*r
    elseif prestoner < 0.998 && stoner > 0.998
        r = 0.5*r
    end
    if r != mixing
        println("we update r. Here r=$r.")
    end
    r
end

# Taken from
# https://discourse.julialang.org/t/non-allocating-matrix-inversion/62264/7
import Base: require_one_based_indexing
import LinearAlgebra: checksquare, chkstride1
import LinearAlgebra.BLAS: @blasfunc, BlasInt
import LinearAlgebra.LAPACK: chkargsok, chklapackerror, liblapack, getrf!, getri!
import LinearAlgebra.LAPACK: geev!

#=
for (getrf, getri, elty) in ((:dgetrf_,:dgetri_,:Float64), (:sgetrf_,:sgetri_,:Float32), (:zgetrf_,:zgetri_,:ComplexF64), (:cgetrf_,:cgetri_,:ComplexF32))
    @eval function getrf!(A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
        require_one_based_indexing(A)
        chkstride1(A)
        chkstride1(ipiv)
        m, n = size(A)
        lda  = max(1,stride(A, 2))
        length(ipiv) ≥ min(m,n) || throw(DimensionMismatch())
        info = Ref{BlasInt}()
        ccall((@blasfunc($getrf), liblapack), Cvoid,
              (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
               Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
              m, n, A, lda, ipiv, info)
        chkargsok(info[])
        A, ipiv, info[] #Error code is stored in LU factorization type
    end
    @eval function getri!(A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt}, work::AbstractVector{$elty})
        require_one_based_indexing(A, ipiv)
        chkstride1(A, ipiv)
        n = checksquare(A)
        if n != length(ipiv)
            throw(DimensionMismatch("ipiv has length $(length(ipiv)), but needs $n"))
        end
        lda = max(1,stride(A, 2))
        lwork = BlasInt(length(work))
        lwork ≥ n || throw(DimensionMismatch())
        info  = Ref{BlasInt}()
        ccall((@blasfunc($getri), liblapack), Cvoid,
              (Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
               Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
              n, A, lda, ipiv, work, lwork, info)
        chklapackerror(info[])
        A
    end
end

struct InvWorkspace{T}
    ipiv::Vector{BlasInt}
    work::Vector{T}
end

function InvWorkspace(A::AbstractMatrix{T}) where {T}
    n = checksquare(A)
    ipiv = Vector{BlasInt}(undef, n)
    work = Vector{T}(undef, 64*n)
    return InvWorkspace{T}(ipiv, work)
end

function inv!(A::AbstractMatrix{T}, W::InvWorkspace{T}) where {T}
    getrf!(A, W.ipiv)
    return getri!(A, W.ipiv, W.work)
end
=#

function getrf!(A::AbstractMatrix{ComplexF64}, ipiv::AbstractVector{BlasInt})
    require_one_based_indexing(A)
    chkstride1(A)
    chkstride1(ipiv)
    m, n = size(A)
    lda  = max(1,stride(A, 2))
    length(ipiv) ≥ min(m,n) || throw(DimensionMismatch())
    info = Ref{BlasInt}()
    ccall((@blasfunc(zgetrf_), liblapack), Cvoid,
          (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
           Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
          m, n, A, lda, ipiv, info)
    chkargsok(info[])
    A, ipiv, info[] #Error code is stored in LU factorization type
end
function getri!(A::AbstractMatrix{ComplexF64}, ipiv::AbstractVector{BlasInt}, work::AbstractVector{ComplexF64})
    require_one_based_indexing(A, ipiv)
    chkstride1(A, ipiv)
    n = checksquare(A)
    if n != length(ipiv)
        throw(DimensionMismatch("ipiv has length $(length(ipiv)), but needs $n"))
    end
    lda = max(1,stride(A, 2))
    lwork = BlasInt(length(work))
    lwork ≥ n || throw(DimensionMismatch())
    info  = Ref{BlasInt}()
    ccall((@blasfunc(zgetri_), liblapack), Cvoid,
          (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}),
          n, A, lda, ipiv, work, lwork, info)
    chklapackerror(info[])
    A
end

struct InvWorkspace{ComplexF64}
    ipiv::Vector{BlasInt}
    work::Vector{ComplexF64}
end

function InvWorkspace(A::AbstractMatrix{ComplexF64}) 
    n = checksquare(A)
    ipiv = Vector{BlasInt}(undef, n)
    work = Vector{ComplexF64}(undef, 64*n)
    return InvWorkspace{ComplexF64}(ipiv, work)
end

function inv!(A::AbstractMatrix{ComplexF64}, W::InvWorkspace{ComplexF64}) 
    getrf!(A, W.ipiv)
    # Return不要?
    return getri!(A, W.ipiv, W.work)
end


function zgeev!(A::AbstractMatrix{ComplexF64}, w::AbstractVector{ComplexF64}, work::AbstractVector{ComplexF64}, rwork::AbstractVector{Float64}, lwork::BlasInt)
    require_one_based_indexing(A)
    chkstride1(A)
    chkstride1(w)
    chkstride1(work)
    chkstride1(rwork)
    jobvl = 'N'
    jobvr = 'N'
    n = size(A,1)
    lda  = max(1,stride(A, 2))
    vl = Vector{ComplexF64}(undef,1)
    ldvl::BlasInt = 1
    vr = Vector{ComplexF64}(undef,1)
    ldvr::BlasInt = 1
    #lwork = BlasInt(length(W.work))
    #lwork ≥ n || throw(DimensionMismatch())
    info = Ref{BlasInt}()
    ccall((@blasfunc(zgeev_), liblapack), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
               Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
               Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
               Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
               Ptr{Float64}, Ref{BlasInt}),
               jobvl, jobvr, n, A, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    chkargsok(info[])
    A, w, info[] #Error code is stored in LU factorization type
end

struct EigWorkspace{ComplexF64}
    lwork::Int64
    work::Vector{ComplexF64}
    rwork::Vector{Float64}
end

function EigWorkspace(A::AbstractMatrix{ComplexF64}, lwork::Int64)
    n = size(A,1)
    work = Vector{ComplexF64}(undef, lwork)
    rwork = Vector{Float64}(undef, 2*n)
    return EigWorkspace{ComplexF64}(lwork,work,rwork)
end

function lwork_check(A::AbstractMatrix{ComplexF64}, w::AbstractVector{ComplexF64})
    W = EigWorkspace(A,1)
    zgeev!(A, w, W.work, W.rwork, -1)
    lwork = Int64(real(W.work[1]))
    return lwork
end


function eig!(A::AbstractMatrix{ComplexF64}, w::AbstractVector{ComplexF64}, W::EigWorkspace{ComplexF64})
    #Return不要？
    return zgeev!(A, w, W.work, W.rwork, W.lwork)
end



"""
Force garbage collection in Julia and Python
"""
function gc_collect()
    GC.gc()
end

function typestable(@nospecialize(f), @nospecialize(t); checkonlyany=false)
    v = code_typed(f, t)
    stable = true
    for vi in v
        for (name, ty) in zip(vi[1].slotnames, vi[1].slottypes)
            !(ty isa Type) && continue
            if (checkonlyany && ty === Any) ||
               (!checkonlyany && (!Base.isdispatchelem(ty) || ty == Core.Box))
                stable = false
                println("Type instability is detected! the variable is $(name) ::$ty")
            end
        end
    end
    return stable
end


"""
a_{ij} ≈ sum_d b_{id} * c_{dj}
"""
function decomp(a::AbstractMatrix, eps)
    svd_res = svd(a)
    idx = svd_res.S ./ svd_res.S[1] .> eps
    sqrt_S = sqrt.(svd_res.S[idx])
    a1 = svd_res.U[:, idx] .* reshape(sqrt_S, 1, :)
    a2 = reshape(sqrt_S, :, 1) .* svd_res.Vt[idx, :]
    return a1, transpose(a2)
end