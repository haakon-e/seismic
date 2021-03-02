module Inversion

# Local packages
using ..Grids

# External packages
using LinearAlgebra

# Export methods
export Rxx, Ryy, Rnn, computeUncertaintyMatrix

# Methods
"""
    Rxx(grid::Grid, σ::T, λ::T) where T<:Real

Compute the second-order moment of the physical variable x
given correlation level σ and de
"""
function Rxx(grid::Grid, σ::T, λ::T) where T<:Real
    n = grid.nx*grid.ny; kernel = zeros(n,n)
    for (i, xi) in enumerate(grid.centers)
        for (j, xj) in enumerate(grid.centers)
            kernel[i,j] = σ^2 * exp(-distance(xi,xj) / λ)
        end
    end
    kernel
end

function Ryy(grid::Grid, σ::T) where T<:Real
    σ^2 * I
end

function Rnn(E, rxx, ryy)
    ryy - E * rxx * E'
end

function computeUncertaintyMatrix(E,rxx,rnn)
    P = rxx - rxx * E' * inv(E * rxx * E' + rnn) * E * rxx
    return P
end


end # module