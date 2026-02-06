# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# Lattice utilities (1D for now, extendable to 2D later)

@inline function shifted_site(i::Int, d::Int, L::Int; periodic::Bool)
    j = i + d
    if periodic
        return mod1(j, L)
    end
    return (1 <= j <= L) ? j : nothing
end


"""
    min_image(d::Int, L::Int) -> Int

Calculate the minimum image convention displacement for periodic boundary conditions.

Maps a displacement `d` to the range `[-L/2, L/2)` using periodic boundary conditions.
For even `L`, this uses a symmetric convention where displacements are minimized
by choosing the nearest periodic image.

# Arguments
- `d::Int`: The displacement to be mapped
- `L::Int`: The length of the periodic domain

# Returns
- `Int`: The minimum image displacement in the range `[-L/2, L/2)`
"""
@inline function min_image(d::Int, L::Int)
    # Map displacement to [-L/2, L/2) (for even L, uses symmetric convention)
    return mod(d + L ÷ 2, L) - (L ÷ 2)
end
