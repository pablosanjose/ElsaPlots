# function plot(m::Mesh{2}; resolution = (1024, 1024), kw...)
#     scene = Scene(resolution = resolution)
#     cam3d!(scene)

#     plot = plot!(scene, m)

#     return plot

# end

# function default_theme(scene::SceneLike, ::Type{<: Plot(Mesh)})
#     Theme(
#         wireframe = 0
#         )
# end

# function AbstractPlotting.plot!(plot::Plot(Mesh))
#     mesh = to_value(plot[1])

#     # for i in eachindex(bm.elements.indices)
#     #     e = bm.elements.indices[i]
#     #     sites = [Point3D(s) + Point3D(SVector(0,0,0.1*i)) for s in bm.mesh.sublats[1].sites]
#     #     mesh!(plot, sites, [e[j] for j in 1:3], color = :white, strokewidth = 8)
#     # end
#     verts = Elsa.vertices(mesh)
#     for group in bm.mesh.elements.groupsintra
#         mesh!(plot, Point3D.(verts), [e[j] for e in group, j in 1:3], color = :white, strokewidth = 8)
#         # plot!(plot, bm.mesh)

#         centers = [padright((sites[e[1]]+ sites[e[2]] + sites[e[3]])/3, Val(3)) for e in group]
#         normals = [cross(padright(sites[e[2]] - sites[e[1]], Val(3)),
#                         padright(sites[e[3]] - sites[e[1]], Val(3))) for e in group]
#         segments = [Point3D(c) => Point3D(c + 60 * n) for (c,n) in zip(centers, normals)]
#         linesegments!(plot, segments, color = :blue, strokewidth = 8)
#     end
#     return plot
#  end



function plot(bs::Bandstructure{2}; resolution = (1024, 1024), kw...) where {T}
    scene = Scene(resolution = resolution)
    cam3d!(scene)
    plot = plot!(scene, bs)
    return plot
end

function default_theme(scene::SceneLike, ::Type{<: Plot(Bandstructure)})
    Theme(
        wireframe = 0,
        colorscheme = map(t -> RGBAf0(t...),
            # Set1_9 from ColorSchemes.jl
            # ((0.894,0.102,0.11), (0.216,0.494,0.722), (0.302,0.686,0.29), (0.596,0.306,0.639), (1.0,0.498,0.0), (1.0,1.0,0.2), (0.651,0.337,0.157), (0.969,0.506,0.749), (0.6,0.6,0.6)))
            # sqrt(Set1_9) from ColorSchemes.jl
            # ((0.946, 0.319, 0.332), (0.465, 0.703, 0.85), (0.55, 0.828, 0.539), (0.772, 0.553, 0.799), (1.0, 0.706, 0.0), (1.0, 1.0, 0.447), (0.807, 0.581, 0.396), (0.984, 0.711, 0.865), (0.775, 0.775, 0.775)))
            ((0.973, 0.565, 0.576), (0.682, 0.838, 0.922), (0.742, 0.91, 0.734), (0.879, 0.744, 0.894), (1.0, 0.84, 0.0), (1.0, 1.0, 0.669), (0.898, 0.762, 0.629), (0.992, 0.843, 0.93), (0.88, 0.88, 0.88)))
            # Mathematica default
            #((0.410,0.067,0.031),(0.860,0.400,0.027),(0.940,0.780,0.000),(0.640,0.760,0.900),(0.310,0.370,0.650),(0.600,0.550,0.810),(0.150,0.051,0.100),(0.870,0.530,0.640),(0.720,0.130,0.250)))
        )
end

function AbstractPlotting.plot!(plot::Plot(Bandstructure))
    bs = to_value(plot[1])
    colors = cycle(plot[:colorscheme][])
    for (band, color) in zip(bs.bands, colors)
        vertices = band.mesh.vertices
        simplices = [s[j] for s in band.mesh.simplices, j in 1:3]
        mesh!(plot, vertices, simplices, color = color, strokewidth = 8)
    end

    # sites = Point3D.(bs.mesh.lattice.sublats[1].sites)
    # for i in 1:nsubbands
    # # i=1
    #     group = bs.mesh.elements.groupsintra[i]
    #     mesh!(plot, sites, [e[j] for e in group, j in 1:3], color = colors[i], strokewidth = 8)
    # end

    # nsubbands = Elsa.ngroups(bs.mesh.elements)
    # colors = collect(take(cycle(plot[:colorscheme][]), nsubbands))

    # sites = Point3D.(bs.mesh.lattice.sublats[1].sites)
    # for (i, subband) in enumerate(bs.mesh.elements.groups)
    #     mesh!(plot, sites, [e[j] for e in subband, j in 1:3], color = colors[i], strokewidth = 8)
    # end

    # @show typeof(plot)
    # plot!(plot, bs.mesh.lattice)

    # sites = bs.mesh.lattice.sublats[1].sites
    # for group in bs.mesh.elements.groupsintra
    #     centers = [(sites[e[1]]+ sites[e[2]] + sites[e[3]])/3 for e in group]
    #     normals = [cross(sites[e[2]] - sites[e[1]], sites[e[3]] - sites[e[1]]) for e in group]
    #     segments = [Point3D(c) => Point3D(c + 60 * n) for (c,n) in zip(centers, normals)]
    #     linesegments!(plot, segments, color = :blue, strokewidth = 8)
    # end

    return plot
 end
