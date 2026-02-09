# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Triple density correlators
# ----------------------------

"""
    precompute_n(psi, sites, opname::String)

Precompute ⟨n_i⟩ for all sites and return a vector.
"""
function precompute_n(psi::MPS, sites, opname::String)
    L = length(sites)
    n = zeros(Float64, L)
    for i in 1:L
        n[i] = expect_n(psi, sites, opname, i)
    end
    return n
end

"""
    precompute_nn(psi, sites, opname::String)

Precompute ⟨n_i n_j⟩ for all site pairs and return a symmetric matrix.
"""
function precompute_nn(psi::MPS, sites, opname::String)
    L = length(sites)
    nn = zeros(Float64, L, L)
    for i in 1:L
        nn[i, i] = expect_nn(psi, sites, opname, i, i)
        for j in (i + 1):L
            v = expect_nn(psi, sites, opname, i, j)
            nn[i, j] = v
            nn[j, i] = v
        end
    end
    return nn
end

@inline function _nn_no(nvec::AbstractVector, nnmat::AbstractMatrix, i::Int, j::Int)
    return i == j ? (nnmat[i, i] - nvec[i]) : nnmat[i, j]
end

"""
    connected_nnn(psi, sites, opname::String, i, j, k)

Connected 3-point density correlator (third cumulant) built from plain moments:

⟨n_i n_j n_k⟩_c =
  ⟨n_i n_j n_k⟩
  - ⟨n_i n_j⟩⟨n_k⟩ - ⟨n_i n_k⟩⟨n_j⟩ - ⟨n_j n_k⟩⟨n_i⟩
  + 2 ⟨n_i⟩⟨n_j⟩⟨n_k⟩
"""
function connected_nnn(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int)
    ni = expect_n(psi, sites, opname, i)
    nj = expect_n(psi, sites, opname, j)
    nk = expect_n(psi, sites, opname, k)

    nij = expect_nn(psi, sites, opname, i, j)
    nik = expect_nn(psi, sites, opname, i, k)
    njk = expect_nn(psi, sites, opname, j, k)

    nijk = expect_nnn(psi, sites, opname, i, j, k)
    return nijk - nij * nk - nik * nj - njk * ni + 2 * ni * nj * nk
end

"""
    expect_nn_no(psi, sites, opname::String, i, j)

Normal-ordered 2-point density correlator:

- if i != j:  ⟨:n_i n_j:⟩ = ⟨n_i n_j⟩
- if i == j:  ⟨:n_i^2:⟩ = ⟨n_i(n_i-1)⟩
"""
function expect_nn_no(psi::MPS, sites, opname::String, i::Int, j::Int)
    if i != j
        return expect_nn(psi, sites, opname, i, j)
    end
    n2 = expect_product(psi, sites, [(opname, i), (opname, i)]) # ⟨n^2⟩
    n1 = expect_n(psi, sites, opname, i)                       # ⟨n⟩
    return n2 - n1                                              # ⟨n(n-1)⟩
end

"""
    expect_nnn_no(psi, sites, opname::String, i, j, k)

Normal-ordered 3-point density correlator:

- all distinct: ⟨:n_i n_j n_k:⟩ = ⟨n_i n_j n_k⟩
- one repeated (e.g. i=j): ⟨:n_i^2 n_k:⟩ = ⟨n_i(n_i-1) n_k⟩
- all equal: ⟨:n_i^3:⟩ = ⟨n_i(n_i-1)(n_i-2)⟩
"""
function expect_nnn_no(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int)
    if i != j && i != k && j != k
        return expect_nnn(psi, sites, opname, i, j, k)
    end

    # All equal
    if i == j && j == k
        n3 = expect_product(psi, sites, [(opname, i), (opname, i), (opname, i)]) # ⟨n^3⟩
        n2 = expect_product(psi, sites, [(opname, i), (opname, i)])              # ⟨n^2⟩
        n1 = expect_n(psi, sites, opname, i)                                     # ⟨n⟩
        return n3 - 3 * n2 + 2 * n1                                               # ⟨n(n-1)(n-2)⟩
    end

    # Exactly one pair equal
    if i == j
        n2k = expect_product(psi, sites, [(opname, i), (opname, i), (opname, k)]) # ⟨n_i^2 n_k⟩
        nk = expect_nn(psi, sites, opname, i, k)                                  # ⟨n_i n_k⟩
        return n2k - nk                                                           # ⟨n_i(n_i-1) n_k⟩
    elseif i == k
        n2j = expect_product(psi, sites, [(opname, i), (opname, i), (opname, j)]) # ⟨n_i^2 n_j⟩
        nj = expect_nn(psi, sites, opname, i, j)                                  # ⟨n_i n_j⟩
        return n2j - nj                                                           # ⟨n_i(n_i-1) n_j⟩
    else
        # j == k
        n2i = expect_product(psi, sites, [(opname, j), (opname, j), (opname, i)]) # ⟨n_j^2 n_i⟩
        ni = expect_nn(psi, sites, opname, i, j)                                  # ⟨n_i n_j⟩
        return n2i - ni                                                           # ⟨n_j(n_j-1) n_i⟩
    end
end

"""
    connected_nnn_no(psi, sites, opname::String, i, j, k)

Connected 3-point correlator built out of *normal-ordered* moments:

C^(3)_{ijk} = ⟨:n_i n_j n_k:⟩
  - ⟨:n_i n_j:⟩⟨n_k⟩ - ⟨:n_i n_k:⟩⟨n_j⟩ - ⟨:n_j n_k:⟩⟨n_i⟩
  + 2 ⟨n_i⟩⟨n_j⟩⟨n_k⟩
"""
function connected_nnn_no(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int)
    ni = expect_n(psi, sites, opname, i)
    nj = expect_n(psi, sites, opname, j)
    nk = expect_n(psi, sites, opname, k)

    nij = expect_nn_no(psi, sites, opname, i, j)
    nik = expect_nn_no(psi, sites, opname, i, k)
    njk = expect_nn_no(psi, sites, opname, j, k)

    nijk = expect_nnn_no(psi, sites, opname, i, j, k)
    return nijk - nij * nk - nik * nj - njk * ni + 2 * ni * nj * nk
end

"""
    expect_nnn_no_cached(psi, sites, opname::String, i, j, k, nvec, nnmat)

Normal-ordered 3-point correlator using cached ⟨n_i⟩ and ⟨n_i n_j⟩ when possible.
"""
function expect_nnn_no_cached(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int,
    nvec::AbstractVector, nnmat::AbstractMatrix)
    if i != j && i != k && j != k
        return expect_nnn(psi, sites, opname, i, j, k)
    end

    # All equal
    if i == j && j == k
        n3 = expect_product(psi, sites, [(opname, i), (opname, i), (opname, i)]) # ⟨n^3⟩
        n2 = nnmat[i, i]                                                        # ⟨n^2⟩
        n1 = nvec[i]                                                            # ⟨n⟩
        return n3 - 3 * n2 + 2 * n1                                              # ⟨n(n-1)(n-2)⟩
    end

    # Exactly one pair equal
    if i == j
        n2k = expect_product(psi, sites, [(opname, i), (opname, i), (opname, k)]) # ⟨n_i^2 n_k⟩
        nk = nnmat[i, k]                                                         # ⟨n_i n_k⟩
        return n2k - nk                                                          # ⟨n_i(n_i-1) n_k⟩
    elseif i == k
        n2j = expect_product(psi, sites, [(opname, i), (opname, i), (opname, j)]) # ⟨n_i^2 n_j⟩
        nj = nnmat[i, j]                                                         # ⟨n_i n_j⟩
        return n2j - nj                                                          # ⟨n_i(n_i-1) n_j⟩
    else
        # j == k
        n2i = expect_product(psi, sites, [(opname, j), (opname, j), (opname, i)]) # ⟨n_j^2 n_i⟩
        ni = nnmat[i, j]                                                         # ⟨n_i n_j⟩
        return n2i - ni                                                          # ⟨n_j(n_j-1) n_i⟩
    end
end

"""
    connected_nnn_no_cached(psi, sites, opname::String, i, j, k, nvec, nnmat)

Connected 3-point correlator built from *normal-ordered* moments using cached
⟨n_i⟩ and ⟨n_i n_j⟩.
"""
function connected_nnn_no_cached(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int,
    nvec::AbstractVector, nnmat::AbstractMatrix)
    ni = nvec[i]
    nj = nvec[j]
    nk = nvec[k]

    nij = _nn_no(nvec, nnmat, i, j)
    nik = _nn_no(nvec, nnmat, i, k)
    njk = _nn_no(nvec, nnmat, j, k)

    nijk = expect_nnn_no_cached(psi, sites, opname, i, j, k, nvec, nnmat)
    return nijk - nij * nk - nik * nj - njk * ni + 2 * ni * nj * nk
end

"""
    transl_avg_nnn(psi, sites, opname::String, r::Int, s::Int; periodic::Bool=false)

Translationally-average the 3-point correlator over the "anchor" site i:

G^(3)(r, s) = (1/N) * Σ_i ⟨n_i n_{i+r} n_{i+s}⟩
"""
function transl_avg_nnn(psi::MPS, sites, opname::String, r::Int, s::Int; periodic::Bool=false)
    L = length(sites)
    acc = 0.0
    n = 0
    for i in 1:L
        j = shifted_site(i, r, L; periodic=periodic)
        j === nothing && continue
        k = shifted_site(i, s, L; periodic=periodic)
        k === nothing && continue
        acc += expect_nnn(psi, sites, opname, i, j, k)
        n += 1
    end
    n == 0 && throw(ArgumentError("no valid anchors for r=$r, s=$s with L=$L (periodic=$periodic)"))
    return acc / n, n
end

"""
    transl_avg_connected_nnn(psi, sites, opname::String, r::Int, s::Int; periodic::Bool=false)

Translationally-average the connected 3-point correlator (third cumulant).
Returns `(C, N)`.
"""
function transl_avg_connected_nnn(psi::MPS, sites, opname::String, r::Int, s::Int; periodic::Bool=false)
    L = length(sites)
    acc = 0.0
    n = 0
    for i in 1:L
        j = shifted_site(i, r, L; periodic=periodic)
        j === nothing && continue
        k = shifted_site(i, s, L; periodic=periodic)
        k === nothing && continue
        acc += connected_nnn(psi, sites, opname, i, j, k)
        n += 1
    end
    n == 0 && throw(ArgumentError("no valid anchors for r=$r, s=$s with L=$L (periodic=$periodic)"))
    return acc / n, n
end

"""
    transl_avg_nnn_no(psi, sites, opname::String, r::Int, s::Int; periodic::Bool=false)

Translationally-average the *normal-ordered* 3-point correlator:

G_no^(3)(r, s) = (1/N) * Σ_i ⟨:n_i n_{i+r} n_{i+s}:⟩
"""
function transl_avg_nnn_no(psi::MPS, sites, opname::String, r::Int, s::Int; periodic::Bool=false)
    L = length(sites)
    acc = 0.0
    n = 0
    for i in 1:L
        j = shifted_site(i, r, L; periodic=periodic)
        j === nothing && continue
        k = shifted_site(i, s, L; periodic=periodic)
        k === nothing && continue
        acc += expect_nnn_no(psi, sites, opname, i, j, k)
        n += 1
    end
    n == 0 && throw(ArgumentError("no valid anchors for r=$r, s=$s with L=$L (periodic=$periodic)"))
    return acc / n, n
end

"""
    transl_avg_connected_nnn_no(psi, sites, opname::String, r::Int, s::Int; periodic::Bool=false)

Translationally-average the connected 3-point correlator using *normal-ordered* moments.
Returns `(C, N)`.
"""
function transl_avg_connected_nnn_no(psi::MPS, sites, opname::String, r::Int, s::Int; periodic::Bool=false)
    L = length(sites)
    acc = 0.0
    n = 0
    for i in 1:L
        j = shifted_site(i, r, L; periodic=periodic)
        j === nothing && continue
        k = shifted_site(i, s, L; periodic=periodic)
        k === nothing && continue
        acc += connected_nnn_no(psi, sites, opname, i, j, k)
        n += 1
    end
    n == 0 && throw(ArgumentError("no valid anchors for r=$r, s=$s with L=$L (periodic=$periodic)"))
    return acc / n, n
end

"""
    transl_avg_connected_nnn_no_cached(psi, sites, opname::String, r::Int, s::Int;
                                       periodic::Bool=false, nvec, nnmat)

Translationally-average the connected 3-point correlator using cached ⟨n_i⟩ and
⟨n_i n_j⟩. Returns `(C, N)`.
"""
function transl_avg_connected_nnn_no_cached(psi::MPS, sites, opname::String, r::Int, s::Int;
    periodic::Bool=false, nvec::AbstractVector, nnmat::AbstractMatrix)
    L = length(sites)
    acc = 0.0
    n = 0
    for i in 1:L
        j = shifted_site(i, r, L; periodic=periodic)
        j === nothing && continue
        k = shifted_site(i, s, L; periodic=periodic)
        k === nothing && continue
        acc += connected_nnn_no_cached(psi, sites, opname, i, j, k, nvec, nnmat)
        n += 1
    end
    n == 0 && throw(ArgumentError("no valid anchors for r=$r, s=$s with L=$L (periodic=$periodic)"))
    return acc / n, n
end
