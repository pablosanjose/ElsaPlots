module ElsaPlots

using Elsa, Makie, GeometryTypes, LinearAlgebra, SparseArrays
import AbstractPlotting: default_theme, Plot, plot!, plot, to_value
import ColorTypes: RGBA
using Base.Iterators: take, cycle
using Statistics: mean
using Elsa: vertices, Hamiltonian, Lattice, Bandstructure, padright, siterange

export plot, scale!

include("tools.jl")
include("plotlattice.jl")
include("plotbands.jl")

end
