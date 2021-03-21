module Visualization

using Plots

using ..Grids

export plot_grid_rays

"""
    plot_grid_rays(grid, src, recs; heatmap_col=:tab20)

Visualize the grid, sources, receivers, and ray paths.

For debugging purposes.
"""
function plot_grid_rays(grid::Grid{I,T}, src::Point{T}, recs::Array{Point{T},1}; heatmap_col=:tab20) where {I, T}
    # Plot grid
    p = heatmap(
        grid.xticks,# .- grid.Δx/2,
        grid.yticks,#.- grid.Δy/2, 
        reshape(1:grid.nx*grid.ny,(grid.nx,grid.ny)), 
        color=heatmap_col,
        cbar=false  # hide colorbar
    );
    
    # Annotate grid numbering
    annts = [(pt.x, pt.y, text(i, 7, :white)) for (i, pt) in enumerate(grid.centers)];
    annotate!(annts...);
    
    # Display receivers (floats)
    rx = [r.x for r in recs]; ry = [r.y for r in recs]
    scatter!(rx,ry,label="floats");
    floats_annot = [(pt.x, pt.y, text(i, 10, :bottom)) for (i, pt) in enumerate(recs)];
    annotate!(floats_annot...);
    
    # Display source
    scatter!([src.x],[src.y],label="source");
    
    # Display rays between source and receivers
    for i in 1:length(recs)
        plot!([src.x, rx[i]], [src.y, ry[i]], color=3, lab=false);
    end

    display(p)
end

end # module