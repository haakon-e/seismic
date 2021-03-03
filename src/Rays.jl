module Rays

using ..Grids

export decompose_ray
"""
For each receiver, find the ray segments (and ids) through the cells
"""
function decompose_ray(
        rec::Point{T},
        src::Point{T},
        grid::G,
    ) where {T<:Real, G}
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

    # Remove any points not within the domain
    filter!(pt -> (
            minimum(grid.xticks) <= pt.x <= maximum(grid.xticks)) && (
            minimum(grid.yticks) <= pt.y <= maximum(grid.yticks)),
        pts)

    # Sort intersection points by distance to source
    # (first point closest to source, last point closest to receiver)
    sort!(pts, by=(p->distance(p, src)))

    # Add `src` to the end of the point list
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
        found = false
        for (ind, node) in pairs(grid.nodes)
            cond = (node.x <= p.x) & (p.x < node.x + grid.Δx) & (node.y <= p.y) & (p.y < node.y + grid.Δy)
            if cond
                found = true
                push!(line_ind, ind)
                break
            end
        end
    end
    return lengths, line_ind
end


# Design matrix
export distanceMatrix
"""
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

end  # module
