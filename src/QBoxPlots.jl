module QBoxPlots

using QBox, Makie, GeometryTypes, LinearAlgebra
import AbstractPlotting: default_theme, Plot, plot!, to_value
using Base.Iterators: take, cycle
using QBox: nsublats, bravaismatrix, boundingboxlat, padright
using SparseArrays: nonzeros

export plot, scale!

include("plotlattice.jl")
include("plotbands.jl")

end
