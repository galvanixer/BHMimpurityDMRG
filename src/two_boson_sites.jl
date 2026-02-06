# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# 1) Custom two-boson site type
# ----------------------------

"""
Create site indices for a 1D chain of length L with a two-species boson Hilbert space
|na, nb> where 0 <= na <= nmax_a and 0 <= nb <= nmax_b.

If conserve_qns=true, each basis state carries particle-count "labels" (quantum numbers)
for species a and b. ITensors then only allows tensor contractions that preserve those
counts, which speeds things up and avoids unphysical mixing. If false, the index is a
plain dense space with no conservation rules.
"""
function two_boson_siteinds(L::Int; nmax_a::Int, nmax_b::Int, conserve_qns::Bool=true)
    return [siteind("TwoBoson"; nmax_a=nmax_a, nmax_b=nmax_b, conserve_qns=conserve_qns) for _ in 1:L]
end


"""
    ITensors.siteind(::SiteType"TwoBoson"; nmax_a::Int, nmax_b::Int, conserve_qns::Bool=true, tags::String="Site")

Constructs an `Index` for a two-boson site with specified maximum occupation numbers for each boson type.

# Arguments
- `nmax_a::Int`: Maximum occupation number for boson type `a`.
- `nmax_b::Int`: Maximum occupation number for boson type `b`.
- `conserve_qns::Bool=true`: Whether to conserve quantum numbers (`QN`). If `true`, the index will be block-sparse with quantum number labels; if `false`, a dense index is returned.
- `tags::String="Site"`: Tags to attach to the index for identification.

# Returns
- An `Index` object representing the local Hilbert space of a two-boson site, optionally with quantum number conservation.

# Details
- If `conserve_qns` is `false`, returns a dense `Index` of dimension `(nmax_a + 1) * (nmax_b + 1)`.
- If `conserve_qns` is `true`, returns a block-sparse `Index` where each basis state `|na, nb⟩` is assigned quantum numbers `("Na", na)` and `("Nb", nb)`.
- The `tags` argument is augmented with information about the site type and maximum occupation numbers.

"""
function ITensors.siteind(::SiteType"TwoBoson"; nmax_a::Int, nmax_b::Int,
    conserve_qns::Bool=true, tags::String="Site")
    d = (nmax_a + 1) * (nmax_b + 1) # total dimension

    base = TagSet(tags)  # tags is a String like "Site"

    fulltags = addtags(base, "TwoBoson")
    fulltags = addtags(fulltags, "nmaxa=$(nmax_a)")
    fulltags = addtags(fulltags, "nmaxb=$(nmax_b)")

    # If not conserving QNs, just return a plain dense Index
    if !conserve_qns
        return Index(d; tags=fulltags)
    end

    # QN blocks: each basis state |na,nb> has QN(("Na",na),("Nb",nb))
    # Dimension of each block = 1
    qnblocks = Vector{Pair{QN,Int}}(undef, d)
    k = 1
    for na in 0:nmax_a
        for nb in 0:nmax_b
            qnblocks[k] = (QN(("Na", na), ("Nb", nb)) => 1)
            k += 1
        end
    end
    return Index(qnblocks...; tags=fulltags)
end

function _get_cutoffs(s::Index)
    ts = String.(tags(s))
    nmax_a = parse(Int, split(only(filter(t -> startswith(t, "nmaxa="), ts)), "=")[2])
    nmax_b = parse(Int, split(only(filter(t -> startswith(t, "nmaxb="), ts)), "=")[2])
    return nmax_a, nmax_b
end

# Utility: map (na,nb) to linear basis index (1-based)
@inline function _state_index(na::Int, nb::Int, nmax_b::Int)
    # Ordering: na major, nb minor
    return na * (nmax_b + 1) + nb + 1
end

function _op_diag(s::Index, f::Function)
    nmax_a, nmax_b = _get_cutoffs(s)
    O = ITensor(prime(s), dag(s))
    for na in 0:nmax_a, nb in 0:nmax_b
        k = _state_index(na, nb, nmax_b)
        O[prime(s) => k, dag(s) => k] = f(na, nb)
    end
    return O
end

function _op_shift(s::Index, delta_a::Int, delta_b::Int, amp::Function)
    nmax_a, nmax_b = _get_cutoffs(s)
    O = ITensor(prime(s), dag(s))
    for na in 0:nmax_a, nb in 0:nmax_b
        na2 = na + delta_a
        nb2 = nb + delta_b
        if 0 <= na2 <= nmax_a && 0 <= nb2 <= nmax_b
            kin = _state_index(na, nb, nmax_b)
            kout = _state_index(na2, nb2, nmax_b)
            O[prime(s) => kout, dag(s) => kin] = amp(na, nb)
        end
    end
    return O
end

function ITensors.op(::OpName"Na", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_diag(s, (na, nb) -> na)
end

function ITensors.op(::OpName"Nb", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_diag(s, (na, nb) -> nb)
end

function ITensors.op(::OpName"A", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_shift(s, -1, 0, (na, nb) -> sqrt(na))
end

function ITensors.op(::OpName"Adag", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_shift(s, 1, 0, (na, nb) -> sqrt(na + 1))
end

function ITensors.op(::OpName"B", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_shift(s, 0, -1, (na, nb) -> sqrt(nb))
end

function ITensors.op(::OpName"Bdag", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_shift(s, 0, 1, (na, nb) -> sqrt(nb + 1))
end

function ITensors.op(::OpName"Id", ::SiteType"TwoBoson", s::Index; kwargs...)
    return _op_diag(s, (na, nb) -> 1.0)
end
