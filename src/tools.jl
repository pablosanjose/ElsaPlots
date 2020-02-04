#######################################################################
# Tools
#######################################################################
function meandist(h::Hamiltonian)
    distsum = 0.0
    num = 0
    ss = Elsa.sites(h.lattice)
    br = h.lattice.bravais.matrix
    for (row, col, dn) in Elsa.eachindex_nz(h)
        if row != col
            num += 1
            rsrc = ss[col]
            rdst = ss[row] + br * dn
            distsum += norm(rsrc - rdst)
        end
    end
    return iszero(num) ? 0.0 : distsum / num
end

function matrixidx(h::AbstractSparseMatrix, row, col)
    for ptr in nzrange(h, col)
        rowvals(h)[ptr] == row && return ptr
    end
    return 0
end

matrixidx(h::DenseMatrix, row, col) = LinearIndices(h)[row, col]

transparent(rgba::RGBAf0, v = 0.5) = RGBAf0(rgba.r, rgba.g, rgba.b, rgba.alpha * v)

function darken(rgba::RGBAf0, v = 0.66)
    r = max(0, min(rgba.r * (1 - v), 1))
    g = max(0, min(rgba.g * (1 - v), 1))
    b = max(0, min(rgba.b * (1 - v), 1))
    RGBAf0(r,g,b,rgba.alpha)
end
function lighten(rgba, v = 0.66)
    darken(rgba, -v)
end