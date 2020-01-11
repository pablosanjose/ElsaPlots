function plot(h::Hamiltonian{<:Lattice}; resolution = (1000, 1000), kw...)
    scene = hamiltonianplot(h; resolution = resolution, kw...)
    plot = scene[end]
    plot[:tooltips][] && addtooltips!(scene, h)
    scale!(scene)
    return scene
end

@recipe(HamiltonianPlot, hamiltonian) do scene
    Theme(
        allintra = false, allcells = true, intralinks = true, interlinks = true,
        shadedsites = false, shadedlinks = true, dimming = 0.75,
        siteradius = 0.12, siteborder = 3, siteborderdarken = 1.0, linkdarken = 0.0,
        linkthickness = 6, linkoffset = 0.99, linkradius = 0.015,
        tooltips = true, digits = 3,
        _tooltips_rowcolhar = Vector{Tuple{Int,Int,Int}}[],
        colorscheme = map(t -> RGBAf0(t...),
            ((0.960,0.600,.327), (0.410,0.067,0.031),(0.940,0.780,0.000),
            (0.640,0.760,0.900),(0.310,0.370,0.650),(0.600,0.550,0.810),
            (0.150,0.051,0.100),(0.870,0.530,0.640),(0.720,0.130,0.250)))
    )
end

function plot!(plot::HamiltonianPlot)
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
        plotsites!(plot, lat, itr, har.dn, n, csrc´)
    end

    # plot links
    for (n, har) in enumerate(h.harmonics)
        iszero(har.dn) || plot[:allcells][] || break
        for (ssrc, csrc) in zip(sublats, colors)
            csrc´ = iszero(har.dn) ? csrc : transparent(csrc, 1 - plot[:dimming][])
            csrc´ = darken(csrc´, plot[:linkdarken][])
            for (sdst, cdst) in zip(sublats, colors)
                itr = Elsa.indicesnonzeros(har, siterange(lat, sdst), siterange(lat, ssrc))
                plotlinks!(plot, lat, itr, har.dn, n, csrc´)
            end
        end
    end

    return plot
end

function plotsites!(plot, lat, srange, dn, n, color)
    allsites = Elsa.sites(lat)
    br = lat.bravais.matrix
    sites = [padright(allsites[i] + br * dn, Val(3)) for i in srange]
    plot[:tooltips][] && (tt = [(site, 0, n) for site in srange])
    if !isempty(sites)
        plot[:shadedsites][] ? plotsites_hi!(plot, sites, color) : plotsites_lo!(plot, sites, color)
        plot[:tooltips][] && push!(plot[:_tooltips_rowcolhar][], tt)
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

function plotlinks!(plot, lat, itr, dn, n, color)
    links = Pair{SVector{3,Float32},SVector{3,Float32}}[]
    plot[:tooltips][] && (tt = Tuple{Int,Int,Int}[])
    sites = Elsa.sites(lat)
    br = lat.bravais.matrix
    for (row, col) in itr
        iszero(dn) && row == col && continue
        rdst = padright(sites[row] + br * dn, Val(3))
        rsrc = padright(sites[col], Val(3))
        rdst = iszero(dn) ? (rdst + rsrc) / 2 : rdst
        rsrc = rsrc + plot[:siteradius][] * plot[:linkoffset][] * normalize(rdst - rsrc)
        push!(links, rsrc => rdst)
        plot[:tooltips][] && push!(tt, (row, col, n))
    end
    if !isempty(links)
        plot[:shadedlinks][] ? plotlinks_hi!(plot, links, color) : plotlinks_lo!(plot, links, color)
        plot[:tooltips][] && push!(plot[:_tooltips_rowcolhar][], tt)
    end
    return plot
end

function plotlinks_lo!(plot, links, color)
    # linesegments!(plot, links; color = darken(color, plot[:siteborderdarken][]), linewidth = 2+plot[:linkthickness][])
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
        Vec3f0((mp .+ 5 .+ (0, 0))..., 0)
    end
    popup = poly!(campixel(scene), poprect, raw = true, color = RGBAf0(1,1,1,0), visible = visible)
    rect = popup[end]
    translate!(popup, Vec3f0(0, 0, 10000))
    text!(popup, " ", textsize = 30, position = textpos, align = (:center, :center),
        color = :black, strokewidth = 4, strokecolor = :white, raw = true, visible = visible)
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
fix_linesegments_bug(idx, subplot::LineSegments) = Int(idx/2) # catches odd idx
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