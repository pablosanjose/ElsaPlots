function plot(sys::System; resolution = (1024, 1024), kw...)
    scene = Scene(resolution = resolution)
    cam = cam3d!(scene)

    plot = plot!(scene, sys)
    scale!(scene)

    # cam = Makie.cameracontrols(plot)
    b1, b2 = Elsa.boundingbox(sys)
    cam.lookat[] = Vec3D((b1 + b2)/2)
    cam.eyeposition[] = Vec3D((b1 + b2)/2) + Vec3f0(0.,0.01,2.) * Float32(normxy(b1 - b2))
    cam.upvector[] = (0.0, 1.0, 0.0)
    update_cam!(plot, cam)

    return plot
end

function default_theme(scene::SceneLike, ::Type{<:Plot(System)})
    Theme(
        allintra = false, allcells = true, intralinks = true, interlinks = true,
        shaded = false, dimming = 0.75,
        meandr = 1.0, siteradius = 0.12, siteborder = 3, siteborderdarken = 1.0,
        linkthickness = 4, linkoffset = 0.99, linkradius = 0.015,
        colorscheme = map(t -> RGBAf0(t...), ((0.410,0.067,0.031),(0.860,0.400,0.027),(0.940,0.780,0.000),(0.640,0.760,0.900),(0.310,0.370,0.650),(0.600,0.550,0.810),(0.150,0.051,0.100),(0.870,0.530,0.640),(0.720,0.130,0.250)))
        )
end

function AbstractPlotting.plot!(plot::Plot(System))
    sys = to_value(plot[1])
    colors = collect(take(cycle(plot[:colorscheme][]), nsublats(sys)))

    bravais = bravaismatrix(sys)
    intrablock = sys.hamiltonian.intra
    celldist0 = bravais * intrablock.ndist

    for block in sys.hamiltonian.inters
        celldist = bravais * block.ndist
        plot[:allintra][] && 
            plotlinks!(plot, rdrs, intrablock, celldist, colors; dimming = plot[:dimming][])
        plot[:interlinks][] && 
            plotlinks!(plot, rdrs, block, celldist0, colors; dimming = plot[:dimming][])
        plot[:allcells][] && 
            plotsites!(plot, rdrs, celldist, colors; dimming = plot[:dimming][])
    end
    plot[:intralinks][] && 
        plotlinks!(plot, sys, intrablock, celldist0, colors; dimming = 0.0)
    plotsites!(plot, sys, celldist0, colors; dimming = 0.0)

    return plot
 end


 function plotsites!(plot, sys, celldist, colors; dimming = 0.0)
    for (sublat, color) in zip(sys.lattice.sublats, colors)
        colordimmed = transparent(color, 1 - dimming)
        sites = [Point3D(celldist + site) for site in sublat.sites]
        plot[:shaded][] ? drawsites_hi!(plot, sites, colordimmed) : drawsites_lo!(plot, sites, colordimmed)
    end
end

function plotlinks!(plot, sys::System{E,L,T,Tv}, block, celldist, colors; dimming = 0.0) where {E,L,T,Tv}
    rdrs = Elsa.uniquelinks(block, sys)
    normlast = rdr -> norm(last(rdr))
    meandr = norm(zero(T))
    for c in CartesianIndices(rdrs)
        rdr = rdrs[c]
        meandr = max(meandr, isempty(rdr) ? meandr : mean(normlast, rdr))
        (s1, s2) = Tuple(c)
        col1, col2 = darken(colors[s1], 0.1), darken(colors[s2], 0.1)
        col1 = transparent(col1, 1 - dimming)
        iszero(celldist) || (col2 = transparent(col2, 1 - dimming))
        plot[:shaded][] ?
            drawlinks_hi!(plot, rdr, celldist, (col1, col2)) :
            drawlinks_lo!(plot, rdr, celldist, (col1, col2))
    end
    
    iszero(meandr) || (plot[:meandr][] = meandr)
    return nothing
end

function drawsites_lo!(plot, sites, color)
    isempty(sites) || scatter!(plot, sites,
        markersize = 2 * plot[:siteradius][] * plot[:meandr][], color = color, 
        strokewidth = plot[:siteborder][],  strokecolor = darken(color, plot[:siteborderdarken][]))
    return nothing
end

function drawsites_hi!(plot, sites, color)
    isempty(sites) || meshscatter!(plot, sites, markersize = plot[:siteradius], color = color)
    return nothing
end

function drawlinks_lo!(plot, rdr, celldist, (col1, col2))
    isempty(rdr) && return nothing
    segments = [fullsegment(celldist + r, dr, plot[:siteradius][] * plot[:linkoffset][]) for (r, dr) in rdr]
    colsegments = collect(take(cycle((col2, col1)), 2 * length(segments)))
    linesegments!(plot, segments, linewidth = plot[:linkthickness][], color = colsegments)
    return nothing
end

function drawlinks_hi!(plot, rdr, celldist, (col1, col2))
    isempty(rdr) && return nothing
    positions = [Point3D(celldist + r) for (r, _) in rdr]
    rotvectors = [Vec3D(dr) for (r, dr) in rdr]
    scales = [Vec3f0(plot[:linkradius][], plot[:linkradius][], norm(dr)/2) for (r, dr) in rdr]
    cylinder = GLNormalMesh(GeometryTypes.Cylinder{3, Float32}(Point3f0(0., 0., 0.), Point3f0(0., 0, 1.0), Float32(1)), 12)
    meshscatter!(plot, positions, marker = cylinder, markersize = scales, rotations = rotvectors, color = col1)
    cylinder = GLNormalMesh(GeometryTypes.Cylinder{3, Float32}(Point3f0(0., 0., 0.), Point3f0(0., 0, -1.0), Float32(1)), 12)
    meshscatter!(plot, positions,  marker = cylinder, markersize = scales, rotations = rotvectors, color = col2)
    return nothing
end

function fullsegment(r, dr, rad)
    dr2 = dr*(1 - 2rad/norm(dr))/2
    return Point3D(r - dr2) => Point3D(r + dr2) # + Point3D(SVector(0,0,0.2rand()))
end

function halfsegment(r, dr, rad)
    dr2 = dr*(1 - 2rad/norm(dr))/2
    return  Vec3D(dr2)
end


Point3D(r::SVector{3,T}) where T = Point3f0(r)
Point3D(r::SVector{N,T}) where {N,T} = Point3f0(padright(r, zero(Float32), Val(3)))
Vec3D(r::SVector{3,T}) where T = Vec3f0(r)
Vec3D(r::SVector{N,T}) where {N,T} = Vec3f0(padright(r, zero(Float32), Val(3)))

normxy(sv::SVector{3}) = norm(sv[1:2])
normxy(sv) = norm(sv)

function darken(rgba::T, v = 0.66) where T
    r = max(0, min(rgba.r * (1 - v), 1))
    g = max(0, min(rgba.g * (1 - v), 1))
    b = max(0, min(rgba.b * (1 - v), 1))
    T(r,g,b,rgba.alpha)
end
function lighten(rgba, v = 0.66)
    darken(rgba, -v)
end
transparent(rgba::T, v = 0.5) where T = T(rgba.r, rgba.g, rgba.b, rgba.alpha * v)
