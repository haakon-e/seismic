module Rays

using ..Grids

export decompose_ray
export distanceMatrix
export receiverPairs
export relativeDistanceMatrix


"""
    decompose_ray(
        rec::Point{T},
        src::Point{T},
        grid::G,
    ) where {T<:Real, G<:AbstractGrid}

For each receiver, find the ray segments (and ids) through the cells.

Returns the tuple `(lengths, line_ind)` where
`lengths` is an array of distances, each entry corresponing to the
indices given by `line_ind`, which refers to indices of `grid` cells.
"""
function decompose_ray(
        rec::Point{T},
        src::Point{T},
        grid::AbstractGrid,
    ) where T<:Real
    # Find coordinates of interesection with grid
    xs = sort([src.x, rec.x])
    ys = sort([src.y, rec.y])
    xind = @. xs[1] < grid.xticks < xs[2]
    yind = @. ys[1] < grid.yticks < ys[2]
    cx = grid.xticks[xind]  # x ticks that ray passes through
    cy = grid.yticks[yind]  # y ticks that ray passes through

    # Compute angle of ray
    θ = atan(rec.y - src.y, rec.x - src.x)

    # Compute points of intersection with gridlines
    cx_y = @. (cx - src.x) * tan(θ) + src.y  # y-coords of ray intersection with vertical gridlines
    cy_x = @. (cy - src.y) / tan(θ) + src.x  # x-coords of ray intersection with horizontal gridlines

    # Turn intersection coordinates into `Point`s
    x_pts = vcat(cx,cy_x)
    y_pts = vcat(cx_y,cy)
    pts = Point{T}[]
    for (x,y) in zip(x_pts,y_pts)
        push!(pts,Point(x,y))
    end

    # Add source to `pts`. If it is outside the domain, `filter` (below) will remove it again.
    push!(pts, src)

    # Remove any points not within the domain
    filter!(pt -> (
            minimum(grid.xticks) <= pt.x <= maximum(grid.xticks)) && (
            minimum(grid.yticks) <= pt.y <= maximum(grid.yticks)),
        pts)

    # Sort intersection points by distance to source
    # (first point closest to source, last point closest to receiver)
    sort!(pts, by=(p->distance(p, src)))

    # Add `rec` to the end of the point list
    push!(pts, rec)

    # Compute lengths and center points of segments
    lengths = T[]
    centers = Point{T}[]
    for i in 2:length(pts)
        push!(lengths, distance(pts[i-1], pts[i]))
        push!(centers, mean_point(pts[i-1], pts[i]))
    end

    # For each line, find which Grid cell it belongs to
    line_ind = Int64[]
    for p in centers
        for (ind, node) in pairs(grid.nodes)
            cond = (node.x <= p.x) & (p.x < node.x + grid.Δx) & (node.y <= p.y) & (p.y < node.y + grid.Δy)
            if cond
                push!(line_ind, ind)
                break
            end
        end
    end
    return lengths, line_ind
end


# Design matrix
"""
    distanceMatrix(
        receivers::Array{Point{T},1},
        src::Point{T},
        grid::Grid,
    ) where {T<:Real}

Construct matrix which describes the distance that each float (rows)
passes through the grid boxes (columns).
"""
function distanceMatrix(
        receivers::Array{Point{T},1},
        src::Point{T},
        grid::Grid,
	) where {T<:Real}
    D = zeros((length(receivers), grid.nx*grid.ny))
    for (i, rec) in pairs(receivers)  # iterate over receivers
        (lengths, inds) = decompose_ray(rec, src, grid)
        D[i, inds] = lengths
    end
    D
end

"""
    receiverPairs(
        receivers::Array{Point{T},1},
        src::Point{T};
        max_degrees = 10,
        max_distance = 500e3
    ) where {T<:Real}

Construct a matrix that maps from receivers (rows)
to receiver pairs (columns). 

Receiver pairs indicate relative measurements, so for each
row (pair), exactly one column is +1 and one column is -1.
"""
function receiverPairs(
        receivers::Array{Point{T},1},
        src::Point{T};
        max_degrees = 10.,
        max_distance = 500e3
    ) where {T<:Real}
    n_recs = length(receivers)
    M = Array{Int}(undef, 0, n_recs)
    for i in 1:n_recs-1
        rec_i = receivers[i]
        for j in i+1:n_recs
            rec_j = receivers[j]
            # Don't consider receiver pairs far apart or with too large angle relative to the source
            if ∠(rec_i, src, rec_j, in_degrees=true) > max_degrees continue end
            if distance(rec_i, rec_j) > max_distance continue end
            
            m = zeros((1, n_recs))
            m[[i, j]] = [1, -1]
            M = [M; m]  # add to mapping
            println("Added floats $i and $j.")
        end
    end
    return M
end

"""
    relativeDistanceMatrix(
        M::T, D::T
    ) where {T<:AbstractMatrix} 

Compute the relative distance matrix E = M * D, where 
M is the receivers-pairs mapping given by `receiverPairs`, and
D is the distance matrix given by `distanceMatrix`.
"""
function relativeDistanceMatrix(
        M::T, D::T
    ) where {T<:AbstractMatrix} 
    M*D 
end

"""
    relativeDistanceMatrix(
        receivers::Array{Point{T},1},
        src::Point{T},
        grid::Grid;
        min_degrees = 10,
        min_distance = 500e3
    ) where {T<:Real}

Construct matrix which describes the relative distance that pairs
of floats (rows) pass through.

each float (rows)
passes through the grid boxes (columns).
"""
function relativeDistanceMatrix(
        receivers::Array{Point{T},1},
        src::Point{T},
        grid::Grid;
        max_degrees = 10.,
        max_distance = 500e3
    ) where {T<:Real}
    M = receiverPairs(receivers, src, max_degrees=max_degrees, max_distance=max_distance)
    D = distanceMatrix(receivers, src, grid)
    return relativeDistanceMatrix(M, D)
end

end  # module
