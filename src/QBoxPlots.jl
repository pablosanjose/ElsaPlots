module QBoxPlots

using Makie, GeometryTypes
import AbstractPlotting: default_theme, Plot, plot!, to_value
using Base.Iterators: take, cycle

export plot

include("plot.jl")

end
