#TO DO Hami!の削除
# ハミルトニアンの定義からmuを除いた。（2分探索等をしやすくするため)
function Hami!(x::Float64, y::Float64, t_pra::Float64, alpha::Float64, H::Array{ComplexF64,2})
    H[1,1] = -2.0*(cos(x)+cos(y)) + 4.0*t_pra*cos(x)*cos(y) 
    H[2,1] = alpha*(-2.0*sin(y) + 4*t_pra*cos(x)*sin(y) - im*(2.0*sin(x) - 4*t_pra*sin(x)*cos(y)))
    H[1,2] = alpha*(-2.0*sin(y) + 4*t_pra*cos(x)*sin(y) + im*(2.0*sin(x) - 4*t_pra*sin(x)*cos(y)))
    H[2,2] = -2.0*(cos(x)+cos(y)) + 4.0*t_pra*cos(x)*cos(y)
end

# 全ての引数経由で情報を渡す
function Hami(x::Float64, y::Float64, nidof::Int64, t_pra::Float64, alpha::Float64)
    H = Array{Complex{Float64}}(undef, nidof, nidof)
    H[1,1] = -2.0*(cos(x)+cos(y)) + 4.0*t_pra*cos(x)*cos(y) 
    H[2,1] = alpha*(-2.0*sin(y) + 4*t_pra*cos(x)*sin(y) - im*(2.0*sin(x) - 4*t_pra*sin(x)*cos(y)))
    H[1,2] = alpha*(-2.0*sin(y) + 4*t_pra*cos(x)*sin(y) + im*(2.0*sin(x) - 4*t_pra*sin(x)*cos(y)))
    H[2,2] = -2.0*(cos(x)+cos(y)) + 4.0*t_pra*cos(x)*cos(y)
    #Hami!(x, y, t_pra, alpha, H)
    return H
end

# 後方互換性
function Hami_back(x::Float64, y::Float64)
    H = Array{Complex{Float64}}(undef, nidof, nidof)
    Hami!(x, y, t_pra, alpha, H)
    return H
end
