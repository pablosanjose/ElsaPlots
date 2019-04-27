module ElsaPlots

using Elsa, Makie, GeometryTypes, LinearAlgebra, SparseArrays
import AbstractPlotting: default_theme, Plot, plot!, to_value
using Base.Iterators: take, cycle
using Statistics: mean
using Elsa: nsublats, bravaismatrix, padright

export plot, scale!

include("plotlattice.jl")
# include("plotbands.jl")

end
