# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Static structure factor
# ----------------------------

"""
    structure_factor_from_nn(nvec, nnmat; connected::Bool=true, factorial_diagonal::Bool=false)

Compute the static structure factor

`S(k) = (1/L) Σ_{i,j} e^{ik(i-j)} M_{ij}`

where `M = nnmat` (raw) or `M = nnmat - nvec*nvec'` (connected).

If `factorial_diagonal=true`, `nnmat` is interpreted as using same-site
factorial convention `⟨n_i(n_i-1)⟩` on the diagonal, and is converted to
plain moments `⟨n_i^2⟩` before forming `S(k)`.

The returned momentum grid is centered around zero:
- even `L`: `m = -L/2, ..., L/2-1`
- odd `L`: `m = -(L-1)/2, ..., (L-1)/2`
with `k = 2πm/L`.
"""
function structure_factor_from_nn(
    nvec::AbstractVector,
    nnmat::AbstractMatrix;
    connected::Bool=true,
    factorial_diagonal::Bool=false
)
    L = length(nvec)
    size(nnmat, 1) == L && size(nnmat, 2) == L ||
        throw(ArgumentError("nnmat must be LxL with L=length(nvec)"))

    mat = factorial_diagonal ? copy(nnmat) : nnmat
    if factorial_diagonal
        for i in 1:L
            mat[i, i] += nvec[i]
        end
    end
    if connected
        mat = mat .- nvec * transpose(nvec)
    end

    mvals = iseven(L) ? collect(-(L ÷ 2):(L ÷ 2 - 1)) : collect(-(L ÷ 2):(L ÷ 2))
    ks = [2 * pi * m / L for m in mvals]
    sk = zeros(Float64, L)

    for (ik, k) in enumerate(ks)
        acc = 0.0 + 0.0im
        for i in 1:L, j in 1:L
            acc += mat[i, j] * exp(im * k * (i - j))
        end
        sk[ik] = real(acc) / L
    end
    return ks, sk
end
