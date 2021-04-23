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

Typically σ is the prior solution variance at a given
physical location, and the correlation between two
locations is given by an exponential kernel over the
prescribed decorrelation scale λ.
Units:
	- λ [length]
	- σ [time / length]
    - Rxx [(time / length)^2]
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
    Rnn(A::AbstractMatrix, σ::Real, σ_indp::Real)

Compute the prior covariance of the data given correlation level σ and
independent noise σ_indp, with A being a slowness measure and noise
transformation.

Typically, σ [length] is the position uncertainty of the floats.
The correlation level Rnn is then expressed as
`pos. uncertainty (σ)` * `wave slowness (A)` = `travel time uncertainty (Rnn)`
The independent noise σ_indp represents noise in waveform correlation or
independent noise travel-time measurements.
Typically, the matrix A is the product of the slowness `s` and a mapping
from receivers to receiver pairs (see also the method `receiverPairs`).
Units:
	- σ [length]
    - σ_indp [time]
    - A [time / length]
    - Rnn [time^2]
"""
function Rnn(A::AbstractMatrix, σ::Real, σ_indp::Real)
    A * (σ^2*I) * A' + σ_indp^2 * I
end

"""
    Rnn(σ_indp::Real)

Compute the prior data covariance assuming no correlation between data points.
Units:
    - σ_indp [time]
    - Rnn [time^2]
"""
function Rnn(σ_indp::Real)
    σ_indp^2 * I
end

"""
    Ryy(E::AbstractMatrix, rxx::AbstractMatrix, rnn::AbstractMatrix)

Compute the prior covariance of the measurements y given
Rnn, Rxx, and E.
"""
function Ryy(E::AbstractMatrix, rxx::AbstractMatrix, rnn::AbstractMatrix)
    rnn + E * rxx * E'
end

"""
    uncertaintyMatrix(E::AbstractMatrix, rxx::AbstractMatrix, rnn::AbstractMatrix)

Compute the uncertainty matrix for the linear relation
	y = Ex + n,
"""
function uncertaintyMatrix(E::AbstractMatrix, rxx::AbstractMatrix, rnn::AbstractMatrix)
    P = rxx - rxx * E' * inv(E * rxx * E' + rnn) * E * rxx
    return P
end

"""
    computeInversion(
        E::AbstractMatrix, 
        rxx::AbstractMatrix, 
        rnn::AbstractMatrix, 
        y::AbstractArray,
    )

Compute the solution x̃ of the inverse problem.
"""
function computeInversion(
        E::AbstractMatrix, 
        rxx::AbstractMatrix, 
        rnn::AbstractMatrix, 
        y::AbstractArray,
    )
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
