module ElsaPlots

using Elsa, Makie, GeometryTypes, LinearAlgebra
import AbstractPlotting: default_theme, Plot, plot!, to_value
using Base.Iterators: take, cycle
using Elsa: nsublats, bravaismatrix, boundingboxlat, padright

export plot, scale!

include("plotlattice.jl")
include("plotbands.jl")

end
