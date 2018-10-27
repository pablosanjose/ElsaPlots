import AbstractPlotting: default_theme, Plot, plot!, to_value
using Base.Iterators: take, cycle
using QBox: nsublats, bravaismatrix, boundingboxlat, padright

function plot(lat::Lattice; resolution = (1024, 1024), kw...)
    scene = Scene(resolution = resolution)
    cam3d!(scene)

    plot = plot!(scene, lat)

    scale!(scene)

    b1, b2 = boundingboxlat(lat)
    lookat = Vec3D((b1 + b2)/2)
    eye = lookat + Vec3f0(0.,0.01,2.)*normxy(b1 - b2)

    update_cam!(scene, eye, lookat, Vec3f0(0,1,0))
    return plot

end

function default_theme(scene::SceneLike, ::Type{<: Plot(Lattice)})
    Theme(
        allintra = false, allcells = true, intralinks = true, interlinks = true,
        shaded = false, dimming = 0.75, 
        siteradius = 0.05, siteborder = 8, siteborderdarken = 1.0,
        linkthickness = 20, linkoffset = 0.99, linkradius = 0.1,
        colorscheme = map(t -> RGBAf0(t...), ((0.410,0.067,0.031),(0.860,0.400,0.027),(0.940,0.780,0.000),(0.640,0.760,0.900),(0.310,0.370,0.650),(0.600,0.550,0.810),(0.150,0.051,0.100),(0.870,0.530,0.640),(0.720,0.130,0.250)))
        )
end

function AbstractPlotting.plot!(plot::Plot(Lattice))
    lat = to_value(plot[1])
    colors = collect(take(cycle(plot[:colorscheme][]), nsublats(lat)))
    
    celldist0 = bravaismatrix(lat) * lat.links.intralink.ndist
    for ilink in lat.links.interlinks
        celldist = bravaismatrix(lat) * ilink.ndist
        plot[:allintra][] && plotlinks!(plot, lat.links.intralink, celldist; dimming = plot[:dimming][])
        plot[:interlinks][] && plotlinks!(plot, ilink, celldist0, colors; dimming = plot[:dimming][])
        plot[:allcells][] && plotsites!(plot, lat, celldist, colors; dimming = plot[:dimming][])
    end
    plot[:intralinks][] && plotlinks!(plot, lat.links.intralink, celldist0, colors; dimming = 0.0)
    plotsites!(plot, lat, celldist0, colors; dimming = 0.0)

    return plot
 end


 function plotsites!(plot, lat, celldist, colors; dimming = 0.0)
    for (sublat, color) in zip(lat.sublats, colors)
        colordimmed = transparent(color, 1 - dimming)
        sites = [Point3D(celldist + site) for site in sublat.sites]
        plot[:shaded][] ? drawsites_hi!(plot, sites, colordimmed) : drawsites_lo!(plot, sites, colordimmed)
    end
end

function plotlinks!(plot, ilink, celldist, colors; dimming = 0.0)
    for ci in CartesianIndices(ilink.slinks)
        i, j = Tuple(ci)
        col1, col2 = darken(colors[j], 0.1), darken(colors[i], 0.1)
        col2 = transparent(col2, 1 - dimming)
        iszero(celldist) || (col1 = transparent(col1, 1 - dimming))
        slink = ilink.slinks[ci]
        plot[:shaded][] ? 
            drawlinks_hi!(plot, slink.rdr, celldist, (col1, col2)) : 
            drawlinks_lo!(plot, slink.rdr, celldist, (col1, col2))
    end
    return nothing
end

function drawsites_lo!(plot, sites, color)
    isempty(sites) || scatter!(plot, sites, 
        markersize = 2*plot[:siteradius][], color = color, strokewidth = plot[:siteborder][],  strokecolor = darken(color, plot[:siteborderdarken][]))
    return nothing
end

function drawsites_hi!(plot, sites, color)
    isempty(sites) || meshscatter!(plot, sites, markersize = plot[:siteradius], color = color)
    return nothing
end

function drawlinks_lo!(plot, rdr, celldist, (col1, col2))
    isempty(rdr) && return nothing
    segments = [fullsegment(celldist + r, dr, plot[:siteradius][] * plot[:linkoffset][]) for (r, dr) in rdr]
    colsegments = collect(take(cycle((col1, col2)), 2 * length(segments)))
    linesegments!(plot, segments, linewidth = plot[:linkthickness][], color = colsegments)
    return nothing
end

function drawlinks_hi!(plot, rdr, celldist, (col1, col2))
    isempty(rdr) && return nothing
    # positions = view(rdr, :, 1)
    positions = [Point3f0(celldist + r) for (r, _) in rdr]
    segments = [halfsegment(r, dr, 0) for (r, dr) in rdr]
    scales = [Vec3f0(plot[:linkradius][], plot[:linkradius][], norm(dr)) for dr in segments]
    cylinder = GLNormalMesh(Makie.Cylinder{3, Float32}(Point3f0(0., 0., 0.), Point3f0(0., 0, 1.), Float32(1)), 12)
    meshscatter!(plot, positions, marker = cylinder, markersize = scales, rotations = segments, color = col2)
    cylinder = GLNormalMesh(Makie.Cylinder{3, Float32}(Point3f0(0., 0., 0.), Point3f0(0., 0, -1.), Float32(1)), 12)
    meshscatter!(plot, positions,  marker = cylinder, markersize = scales, rotations = segments, color = col1)
    return nothing
end

function fullsegment(r, dr, rad) 
    dr2 = dr*(1 - 2rad/norm(dr))/2
    return Point3D(r - dr2) => Point3D(r + dr2)
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