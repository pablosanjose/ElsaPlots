function plot(h::Hamiltonian{<:Lattice{3}}; kw...)
    scene = hamiltonianplot3d(h; kw...)
    plot = scene[end]
    plot[:tooltips][] && addtooltips!(scene, h)
    # scale!(scene)
    return scene
end

plot(bs::Hamiltonian{<:Lattice{2}}; kw...) = hamiltonianplot2d(bs; kw...)
plot(bs::Hamiltonian{<:Lattice{1}}; kw...) = hamiltonianplot2d(bs; kw...)

plot(bs::Lattice{3}; kw...) = latticeplot3d(bs; kw...)
plot(bs::Lattice{2}; kw...) = latticeplot2d(bs; kw...)
plot(bs::Lattice{1}; kw...) = latticeplot2d(bs; kw...)

@recipe(HamiltonianPlot3D, hamiltonian) do scene
    Theme(
        allintra = false, allcells = true, intralinks = true, interlinks = true,
        shaded = false, dimming = 0.75,
        siteradius = 0.12, siteborder = 3, siteborderdarken = 1.0,
        linkthickness = 4, linkoffset = 0.99, linkradius = 0.015,
        tooltips = true, digits = 3,
        _tooltips_rowcolhar = Vector{Tuple{Int,Int,Int}}[],
        colorscheme = map(t -> RGBAf0(t...),
            ((0.960,0.600,.327), (0.410,0.067,0.031),(0.940,0.780,0.000),
            (0.640,0.760,0.900),(0.310,0.370,0.650),(0.600,0.550,0.810),
            (0.150,0.051,0.100),(0.870,0.530,0.640),(0.720,0.130,0.250)))
    )
end

function plot!(plot::HamiltonianPlot3D)
    h = to_value(plot[1])
    lat = h.lattice
    colors = cycle(plot[:colorscheme][])
    sublats = Elsa.sublats(lat)
    plot[:siteradius][] *= meandist(h)

    # plot sites
    for (n, har) in enumerate(h.harmonics), (ssrc, csrc) in zip(sublats, colors)
        iszero(har.dn) || plot[:allcells][] || break
        csrc´ = iszero(har.dn) ? csrc : transparent(csrc, 1 - plot[:dimming][])
        itr = siterange(lat, ssrc)
        if plot[:tooltips][] 
            t = [(site, 0, n) for site in itr]
            isempty(t) || push!(plot[:_tooltips_rowcolhar][], t)
        end
        plotsites!(plot, lat, itr, har.dn, csrc´)
    end

    # plot links
    for (n, har) in enumerate(h.harmonics)
        iszero(har.dn) || plot[:allcells][] || break
        for (ssrc, csrc) in zip(sublats, colors)
            csrc´ = iszero(har.dn) ? csrc : transparent(csrc, 1 - plot[:dimming][])
            for (sdst, cdst) in zip(sublats, colors)
                itr = Elsa.indicesnonzeros(har, siterange(lat, sdst), siterange(lat, ssrc))
                if plot[:tooltips][]
                    t = [(row, col, n) for (row, col) in itr if !(iszero(har.dn) && row == col)]
                    isempty(t) || push!(plot[:_tooltips_rowcolhar][], t)
                end
                plotlinks!(plot, lat, itr, har.dn, csrc´)
            end
        end
    end

    return plot
end

function plotsites!(plot, lat, srange, dn, color)
    allsites = Elsa.sites(lat)
    br = lat.bravais.matrix
    sites = [allsites[i] + br * dn for i in srange]
    if !isempty(sites)
        plot[:shaded][] ? plotsites_hi!(plot, sites, color) : plotsites_lo!(plot, sites, color)
    end
    return plot
end

function plotsites_lo!(plot, sites, color)
    scatter!(plot, sites;
        color = color,
        markersize = 2 * plot[:siteradius][],
        strokewidth = plot[:siteborder][],
        strokecolor = darken(color, plot[:siteborderdarken][]))
    return nothing
end

function plotsites_hi!(plot, sites, color)
    meshscatter!(plot, sites;
        color = color,
        markersize = plot[:siteradius][] * plot[:linkoffset][])
    return nothing
end

function plotlinks!(plot, lat, itr, dn, color)
    links = Pair{SVector{3,Float32},SVector{3,Float32}}[]
    sites = Elsa.sites(lat)
    br = lat.bravais.matrix
    for (row, col) in itr
        iszero(dn) && row == col && continue
        rsrc = padright(sites[col], Val(3))
        rdst = padright(sites[row] + br * dn, Val(3))
        push!(links, rsrc => iszero(dn) ? (rdst + rsrc) / 2 : rdst)
    end
    if !isempty(links)
        plot[:shaded][] ? plotlinks_hi!(plot, links, color) : plotlinks_lo!(plot, links, color)
    end
    return plot
end

function plotlinks_lo!(plot, links, color)
    linesegments!(plot, links; color = color, linewidth = plot[:linkthickness][])
    return nothing
end

function plotlinks_hi!(plot, links, color)
    positions = [(r1 + r2) / 2 for (r1, r2) in links]
    rotvectors = [r2 - r1 for (r1, r2) in links]
    radius = plot[:linkradius][]
    scales = [Vec3f0(radius, radius, norm(r2 - r1)/2) for (r1, r2) in links]
    cylinder = GeometryTypes.Cylinder(Point3f0(0., 0., -1.0), Point3f0(0., 0, 1.0), Float32(1))
    meshscatter!(plot, positions;
        color = color, marker = cylinder, markersize = scales, rotations = rotvectors)
    return nothing
end

function addtooltips!(scene, h)
    sceneplot = scene[end]
    visible = node(:visible, false)
    N = Elsa.blockdim(h)
    poprect = lift(scene.events.mouseposition) do mp
        FRect((mp .+ 5), 1,1)
    end
    textpos = lift(scene.events.mouseposition) do mp
        Vec3f0((mp .+ 5)..., 0)
    end
    popup = poly!(campixel(scene), poprect, raw = true, color = RGBAf0(1,1,1,0), visible = visible)
    rect = popup[end]
    text!(popup, " ", textsize = 30, position = textpos, color = :black, align = (:center, :center), raw = true, visible = visible)
    text_field = popup[end]
    on(scene.events.mouseposition) do event
        subplot, idx = mouse_selection(scene)
        layer = findfirst(isequal(subplot), sceneplot.plots)
        if layer !== nothing && idx > 0
            idx´ = fix_linesegments_bug(idx, subplot)
            txt = popuptext(sceneplot, layer, idx´, h)
            text_field[1] = txt
            visible[] = true
        else
            visible[] = false
        end
        return
    end
    return scene
end

# idx returned by mouse_selection seems wrong by a factor 2 in LineSegments subplot
fix_linesegments_bug(idx, subplot::LineSegments) = idx ÷ 2
fix_linesegments_bug(idx, subplot) = idx

function popuptext(sceneplot, layer, idx, h)
    cache = sceneplot[:_tooltips_rowcolhar][]
    checkbounds(Bool, cache, layer) || return string("Bug in layer: ", layer, " of ", length(cache))
    checkbounds(Bool, cache[layer], idx) || return string("Bug in idx: ", idx, " of ", length(cache[layer]))
    (row, col_or_zero, haridx) = cache[layer][idx]
    if col_or_zero == 0
        col = iszero(col_or_zero) ? row : col_or_zero
        har = h.harmonics[1]
    else
        col = col_or_zero
        har = h.harmonics[haridx]
    end
    element = round.(har.h[row, col], digits = sceneplot[:digits][])
    isreal = all(o -> imag(o) ≈ 0, element)
    txt = isreal ? matrixstring(real.(element)) : matrixstring(element)
    if col_or_zero == 0
        txt´ = string("Onsite :", txt)
    else
        txt´ = string("Hopping :", txt)
    end
    return txt´
end

matrixstring(x::Number) = string(x)
function matrixstring(s::SMatrix)
    ss = repr("text/plain", s)
    pos = findfirst(isequal('\n'), ss)
    return pos === nothing ? ss : ss[pos:end]
end

# function popuptext(plot, layer, idx, h)
#     nsubs = Elsa.nsublats(h.lattice)
#     totalhar = plot[:allcells][] ? length(h.harmonics) : 1
#     if layer <= nsubs * totalhar
#         txt = popuptext_site(plot, layer, idx, h)
#     else
#         layer´ = layer - nsubs * totalhar
#         txt = popuptext_hop(plot, layer´, idx, h)
#     end
#     return txt
# end

# function popuptext_site(plot, layer, idx, h)
#     nsubs = Elsa.nsublats(h.lattice)
#     sub = mod1(layer, nsubs)
#     siteidx = Elsa.siteindex(h.lattice, sub, idx)
#     nhar = fld(layer - 1, nsubs) + 1
#     har = h.harmonics[nhar]
#     site = Elsa.sites(h.lattice, sub)[idx]
#     onsite = round.(har.h[siteidx, siteidx], digits = plot[:digits][])
#     isreal = all(o -> imag(o) ≈ 0, onsite)
#     txt = isreal ? matrixstring(real.(onsite)) : matrixstring(onsite)
#     txt´ = string("Onsite :\n", txt)
#     return txt´
# end

# function popuptext_hop(plot, layer´, idx, h)
#     nsubs = Elsa.nsublats(h.lattice)
#     totalhar = plot[:allcells][] ? length(h.harmonics) : 1
#     subdst, subsrc, nhar = Tuple(CartesianIndices((1:nsubs, 1:nsubs, 1:totalhar))[layer´])
#     har = h.harmonics[nhar]
#     hopping = round.(nonzeros(har.h)[idx], digits = plot[:digits][])
#     # hopping = round.(har.h[subdst, subsrc], digits = plot[:digits][])
#     isreal = all(o -> imag(o) ≈ 0, hopping)
#     txt = isreal ? matrixstring(real.(hopping)) : matrixstring(hopping)
#     txt´ = string("Hopping :\n", txt)
#     return txt´
# end

# function matrixstring(s::SMatrix{N,M}) where {N,M}
#     ss = string.(Vector(vec(transpose(s))))
#     len = maximum(length, ss) + 1
#     map!(z -> lpad(z, len), ss, ss)
#     vs = vec(ss)
#     counter = 0
#     for row in 1:N, col in 1:M
#         counter += 1
#         col == 1 && (insert!(vs, counter, "["); counter += 1)
#         col == M && row < N && (insert!(vs, counter + 1, "]\n"); counter += 1)
#         col == M && row == N && (insert!(vs, counter + 1, "]"); counter += 1)
#     end
#     return join(vs)
# end

# _cyl(r1, r2, radius) = GeometryTypes.Cylinder(Point3f0(r1), Point3f0(r2), Float32(radius))

# function plot(sys::System; resolution = (1024, 1024), kw...)
#     scene = Scene(resolution = resolution)
#     cam = cam3d!(scene)

#     plot = plot!(scene, sys)
#     scale!(scene)

#     # cam = Makie.cameracontrols(plot)
#     b1, b2 = Elsa.boundingbox(sys)
#     cam.lookat[] = Vec3D((b1 + b2)/2)
#     cam.eyeposition[] = Vec3D((b1 + b2)/2) + Vec3f0(0.,0.01,2.) * Float32(normxy(b1 - b2))
#     cam.upvector[] = (0.0, 1.0, 0.0)
#     update_cam!(plot, cam)

#     return plot
# end

# function default_theme(scene::SceneLike, ::Type{<:Plot(System)})
#     Theme(
#         allintra = false, allcells = true, intralinks = true, interlinks = true,
#         shaded = false, dimming = 0.75,
#         siteradius = 0.12, siteborder = 3, siteborderdarken = 1.0,
#         linkthickness = 4, linkoffset = 0.99, linkradius = 0.015,
#         colorscheme = map(t -> RGBAf0(t...), ((0.960,0.600,.327), (0.410,0.067,0.031),(0.940,0.780,0.000),(0.640,0.760,0.900),(0.310,0.370,0.650),(0.600,0.550,0.810),(0.150,0.051,0.100),(0.870,0.530,0.640),(0.720,0.130,0.250)))
#         )
# end

# function AbstractPlotting.plot!(plot::Plot(System))
#     sys = to_value(plot[1])
#     colors = collect(take(cycle(plot[:colorscheme][]), nsublats(sys)))

#     plot[:siteradius][] *= meandist(sys)

#     bravais = bravaismatrix(sys)
#     intrablock = sys.hamiltonian.intra
#     celldist0 = bravais * intrablock.ndist

#     for block in sys.hamiltonian.inters
#         celldist = bravais * block.ndist
#         plot[:allintra][] &&
#             plotlinks!(plot, sys, intrablock, celldist, colors; dimming = plot[:dimming][])
#         plot[:interlinks][] &&
#             plotlinks!(plot, sys, block, celldist0, colors; dimming = plot[:dimming][])
#         plot[:allcells][] &&
#             plotsites!(plot, sys, celldist, colors; dimming = plot[:dimming][])
#     end

#     plot[:intralinks][] &&
#         plotlinks!(plot, sys, intrablock, celldist0, colors; dimming = 0.0)
#     plotsites!(plot, sys, celldist0, colors; dimming = 0.0)

#     return plot
#  end


#  function plotsites!(plot, sys, celldist, colors; dimming = 0.0)
#     for (sublat, color) in zip(sys.lattice.sublats, colors)
#         colordimmed = transparent(color, 1 - dimming)
#         sites = [Point3D(celldist + site) for site in sublat.sites]
#         plot[:shaded][] ? drawsites_hi!(plot, sites, colordimmed) : drawsites_lo!(plot, sites, colordimmed)
#     end
# end

# function plotlinks!(plot, sys::System{E,L,T,Tv}, block, celldist, colors; dimming = 0.0) where {E,L,T,Tv}
#     rdrs = uniquelinks(block, sys)
#     for c in CartesianIndices(rdrs)
#         rdr = rdrs[c]
#         (s1, s2) = Tuple(c)
#         col1, col2 = darken(colors[s1], 0.1), darken(colors[s2], 0.1)
#         col1 = transparent(col1, 1 - dimming)
#         iszero(celldist) || (col2 = transparent(col2, 1 - dimming))
#         plot[:shaded][] ?
#             drawlinks_hi!(plot, rdr, celldist, (col1, col2)) :
#             drawlinks_lo!(plot, rdr, celldist, (col1, col2))
#     end

#     return nothing
# end

# function drawsites_lo!(plot, sites, color)
#     isempty(sites) || scatter!(plot, sites,
#         markersize = 2 * plot[:siteradius][], color = color,
#         strokewidth = plot[:siteborder][],
#         strokecolor = darken(color, plot[:siteborderdarken][]))
#     return nothing
# end

# function drawsites_hi!(plot, sites, color)
#     isempty(sites) || meshscatter!(plot, sites,
#         markersize = plot[:siteradius][]* plot[:linkoffset][], color = color)
#     return nothing
# end

# function drawlinks_lo!(plot, rdr, celldist, (col1, col2))
#     isempty(rdr) && return nothing
#     segments = [fullsegment(celldist + r, dr, plot[:siteradius][] * plot[:linkoffset][])
#                 for (r, dr) in rdr]
#     colsegments = collect(take(cycle((col2, col1)), 2 * length(segments)))
#     linesegments!(plot, segments, linewidth = plot[:linkthickness][], color = colsegments)
#     return nothing
# end

# function drawlinks_hi!(plot, rdr, celldist, (col1, col2))
#     isempty(rdr) && return nothing
#     positions = [Point3D(celldist + r) for (r, _) in rdr]
#     rotvectors = [Vec3D(dr) for (r, dr) in rdr]
#     scales = [Vec3f0(plot[:linkradius][], plot[:linkradius][], norm(dr)/2) for (r, dr) in rdr]
#     cylinder = GLNormalMesh(GeometryTypes.Cylinder{3, Float32}(Point3f0(0., 0., 0.), Point3f0(0., 0, 1.0), Float32(1)), 12)
#     meshscatter!(plot, positions, marker = cylinder, markersize = scales, rotations = rotvectors, color = col1)
#     cylinder = GLNormalMesh(GeometryTypes.Cylinder{3, Float32}(Point3f0(0., 0., 0.), Point3f0(0., 0, -1.0), Float32(1)), 12)
#     meshscatter!(plot, positions,  marker = cylinder, markersize = scales, rotations = rotvectors, color = col2)
#     return nothing
# end

# function fullsegment(r, dr, rad)
#     dr2 = (dr/2) * (1 - rad / norm(dr/2))
#     return Point3D(r - dr2) => Point3D(r + dr2) # + Point3D(SVector(0,0,0.2rand()))
# end


# Point3D(r::SVector{3,T}) where T = Point3f0(r)
# Point3D(r::SVector{N,T}) where {N,T} = Point3f0(padright(r, zero(Float32), Val(3)))
# Vec3D(r::SVector{3,T}) where T = Vec3f0(r)
# Vec3D(r::SVector{N,T}) where {N,T} = Vec3f0(padright(r, zero(Float32), Val(3)))

# normxy(sv::SVector{3}) = norm(sv[1:2])
# normxy(sv) = norm(sv)

# function darken(rgba::T, v = 0.66) where T
#     r = max(0, min(rgba.r * (1 - v), 1))
#     g = max(0, min(rgba.g * (1 - v), 1))
#     b = max(0, min(rgba.b * (1 - v), 1))
#     T(r,g,b,rgba.alpha)
# end
# function lighten(rgba, v = 0.66)
#     darken(rgba, -v)
# end
