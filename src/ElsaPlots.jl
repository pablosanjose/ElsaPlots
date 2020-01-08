module ElsaPlots

using Elsa, Makie, GeometryTypes, LinearAlgebra, SparseArrays
import AbstractPlotting: default_theme, Plot, plot!, plot, to_value
using Base.Iterators: take, cycle
using Statistics: mean
using Elsa: vertices, Bandstructure

export plot, scale!

# include("interface.jl")
# include("plotlattice.jl")
include("plotbands.jl")

end
