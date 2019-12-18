plot(bs::Bandstructure{1}; kw...) = bandplot2d(bs; kw...)
plot(bs::Bandstructure{2}; kw...) = bandplot3d(bs; kw...)

@recipe(BandPlot2D, bandstructure) do scene
    Theme(
    linewidth = 3,
    colorscheme = map(t -> RGBAf0((0.8 .* t)...),
        ((0.973, 0.565, 0.576), (0.682, 0.838, 0.922), (0.742, 0.91, 0.734),
         (0.879, 0.744, 0.894), (1.0, 0.84, 0.0), (1.0, 1.0, 0.669),
         (0.898, 0.762, 0.629), (0.992, 0.843, 0.93), (0.88, 0.88, 0.88)))
    )
end

function plot!(plot::BandPlot2D)
    bs = to_value(plot[1])
    bands = haskey(plot, :bands) ? to_value(plot[:bands]) : eachindex(bs.bands)
    colors = cycle(plot[:colorscheme][])
    for (nb, color) in zip(bands, colors)
        band = bs.bands[nb]
        vertices = band.mesh.vertices
        lines!(plot, first.(vertices), last.(vertices), linewidth = plot[:linewidth], color = color)
    end
    return plot
 end

@recipe(BandPlot3D, bandstructure) do scene
    Theme(
    wireframe = 0,
    colorscheme = map(t -> RGBAf0(t...),
        ((0.973, 0.565, 0.576), (0.682, 0.838, 0.922), (0.742, 0.91, 0.734),
         (0.879, 0.744, 0.894), (1.0, 0.84, 0.0), (1.0, 1.0, 0.669),
         (0.898, 0.762, 0.629), (0.992, 0.843, 0.93), (0.88, 0.88, 0.88)))
    )
end

function plot!(plot::BandPlot3D)
    bs = to_value(plot[1])
    bands = haskey(plot, :bands) ? to_value(plot[:bands]) : eachindex(bs.bands)
    colors = cycle(plot[:colorscheme][])
    for (nb, color) in zip(bands, colors)
        band = bs.bands[nb]
        vertices = band.mesh.vertices
        simplices = [s[j] for s in band.simplices, j in 1:3]
        if isempty(simplices)
            scatter!(plot, vertices, color = color)
        else
            mesh!(plot, vertices, simplices, color = color, transparency = false)
            linesegments!(plot, (t -> vertices[first(t)] => vertices[last(t)]).(simplices),
                          linewidth = 1)
        end
    end
    return plot
 end