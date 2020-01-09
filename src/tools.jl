#######################################################################
# Tools
#######################################################################
function meandist(h::Hamiltonian)
    distsum = 0.0
    num = 0
    ss = Elsa.sites(h.lattice)
    br = h.lattice.bravais.matrix
    for (row, col, dn) in Elsa.indicesnonzeros(h)
        if row != col
            num += 1
            rsrc = ss[col]
            rdst = ss[row] + br * dn
            distsum += norm(rsrc - rdst)
        end
    end
    return distsum / num
end

transparent(rgba::T, v = 0.5) where T = T(rgba.r, rgba.g, rgba.b, rgba.alpha * v)

function darken(rgba::T, v = 0.66) where T
    r = max(0, min(rgba.r * (1 - v), 1))
    g = max(0, min(rgba.g * (1 - v), 1))
    b = max(0, min(rgba.b * (1 - v), 1))
    T(r,g,b,rgba.alpha)
end
function lighten(rgba, v = 0.66)
    darken(rgba, -v)
end