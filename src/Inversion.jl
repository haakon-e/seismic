"""
Solve the inverse problem where an unknown state vector x is
linearly related to a set of observations y through the relation
	Ex + n = y
where E is the design matrix and n is measurement noise not.
The result is an application of the Gauss-Markov theorem which
produces the values x̃ that minimizes the variance of all possible
linear estimators.

The implementation follows chapter 2.7.2 of
Wunsch, 2009: Discrete Inverse and State Estimation Problems.
Specifically, refer to Eqs. 2.400, 2.403, 2.405.
"""
module Inversion

# Local packages
using ..Grids

# External packages
using LinearAlgebra

# Export methods
export Rxx, Ryy, Rnn 
export uncertaintyMatrix, computeInversion

# Methods
"""
    Rxx(grid::Grid, σ::T, λ::T) where T<:Real

Compute the prior covariance of the physical variable x
given correlation level σ and decorrelation scale λ.

Typically, this is the spatial correlation between
locations in the grid space.
Units:
	- λ [length]
	- σ [time / length]
"""
function Rxx(grid::AbstractGrid, σ::T, λ::T) where T <: Real
    n = grid.nx*grid.ny; kernel = zeros(n,n)
    for (i, xi) in enumerate(grid.centers)
        for (j, xj) in enumerate(grid.centers)
            kernel[i,j] = σ^2 * exp(-distance(xi,xj) / λ)
        end
    end
    kernel
end

"""
    Rnn(A::AbstractMatrix, σ::Real)

Compute the prior covariance of the data given correlation level σ.

Typically, this is the position uncertainty of the floats.
The correlation level σ is typically expressed as
`pos. uncertainty` / `wave speed` = `travel time uncertainty` 
Units:
	- σ [time]
"""
function Rnn(A::AbstractMatrix, σ::Real)
    σ^2 * A * A'
end

"""
    Ryy(E::T, rxx::T, rnn::T) where T <: AbstractMatrix

Compute the prior covariance of the measurements y given
Rnn, Rxx, and E.
"""
function Ryy(E::T, rxx::T, rnn::T) where T <: AbstractMatrix
    rnn + E * rxx * E'
end

"""
    uncertaintyMatrix(E::T, rxx::T, rnn::T) <: AbstractMatrix

Compute the uncertainty matrix for the linear relation
	y = Ex + n,
"""
function uncertaintyMatrix(E::T, rxx::T, rnn::T) where T <: AbstractMatrix
    P = rxx - rxx * E' * inv(E * rxx * E' + rnn) * E * rxx
    return P
end

"""
    computeInversion(E::T, rxx::T, rnn::T, y::AbstractArray) where T<:AbstractMatrix

Compute the solution x̃ of the inverse problem.
"""
function computeInversion(E::T, rxx::T, rnn::T, y::AbstractArray) where T <: AbstractMatrix
	rxx * E' * inv(E * rxx * E' + rnn) * y
end

#
# Unused methods:
#

"""
    noiseProjection(
        receivers::Array{Point{T},1},
        src::Point{T},
    ) where T<:Real

Compute the projection of the position noise onto the
normalized x-y plane.
"""
function noiseProjection(
        receivers::Array{Point{T},1},
        src::Point{T},
    ) where T <: Real
    n_recs = length(receivers)
    B = zeros((n_recs, 2*n_recs))
    for (i, rec) in pairs(receivers)
        θ = atan(rec.y - src.y, rec.x - src.x)
        B[i, [2i-1, 2i]] = [cos(θ), sin(θ)]
    end
    B
end

end # module
