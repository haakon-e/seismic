module Grids

export Point, mean_point, distance
export AbstractGrid, Grid, makeGrid
export bounding_box
export grid_neighbors

"""
A `Point` is a pair of coordinates (x,y).
"""
struct Point{T<:Real}
    x::T
    y::T
end

# Basic operations on pairs of points: Compute the distance between them and their midpoint.
distance(p::Point, q::Point) = sqrt((p.x-q.x)^2 + (p.y-q.y)^2)
mean_point(p::Point, q::Point) = Point((p.x+q.x)/2, (p.y+q.y)/2);


abstract type AbstractGrid end

"""
A `Grid` is a regular 2D grid defined by an array of `nodes`, along with the spacing between the nodes
determined by `Δx` and `Δy`. The nodes define the lower-left corners of the cells of the `Grid`.
Associated to the `nodes` are `centers`, which are the coordinates of the cell centers for each grid point.
For convencience, `xticks` and `yticks` store the x- and y- coordinates in separate arrays, like ticks
on a plot. In addition, `nx` and `ny` is the number of cells along the x- and y-dimension.
"""
mutable struct Grid{I<:Integer,F<:Real} <: AbstractGrid
    nodes::Array{Point{F},1}
    centers::Array{Point{F},1}
    Δx::F
    Δy::F
    xticks::Array{F}
    yticks::Array{F}
    nx::I
    ny::I
end

function makeGrid(nx::I, ny::I, Δx::F, Δy::F) where {I<:Integer, F<:Real}
    xticks = collect(0:Δx:nx*Δx)
    yticks = collect(0:Δy:ny*Δy)
    nodes = Point{Float64}[]
    centers = Point{Float64}[]
    for row in 1:nx
        for col in 1:ny
            x = xticks[row]; y = yticks[col]
            push!(nodes, Point(x, y))
            push!(centers, Point(x + Δx/2, y + Δy/2))
        end
    end
    return Grid(nodes, centers, Δx, Δy, xticks, yticks, nx, ny)
    #new(nodes,centers,Δx,Δy,xticks,yticks,nx,ny)
end


"""
Compute the grid bounding box.
Returns the minimum and maximum x- and y-coordinates of the grid.
"""
function bounding_box(grid::Grid)
    xmin = minimum(grid.xticks); xmax = maximum(grid.xticks)
    ymin = minimum(grid.yticks); ymax = maximum(grid.yticks)
    return [xmin, ymin, xmax, ymax]
end


""""
Returns a matrix where a non-zero value implies that the row (cell) and column (other cell) are neighbors.
"""
function grid_neighbors(grid::Grid)
    c = grid.centers; nc = length(c)
    neighbors = zeros(Int64, nc, nc)
    for (idx, cell) in enumerate(c)
        for odx in idx+1:nc  # Only compare points that haven't been compared before
            ocell = c[odx];
            xngb = (abs(cell.x - ocell.x) ≈ grid.Δx) & (cell.y ≈ ocell.y)
            yngb = (abs(cell.y - ocell.y) ≈ grid.Δy) & (cell.x ≈ ocell.x)
            if xngb | yngb
                neighbors[idx, odx] = neighbors[odx, idx] = 1
            end
        end
    end
    neighbors
end

end  # module
