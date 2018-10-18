module QBoxPlots

using QBox, Makie, GeometryTypes
import AbstractPlotting: default_theme, Plot, plot!, to_value
using Base.Iterators: take, cycle

export plot

include("plotlattice.jl")

end
