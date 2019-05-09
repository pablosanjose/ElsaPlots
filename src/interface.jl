#######################################################################
# Functions to extend the interface with Elsa.jl
#######################################################################
function meandist(sys)
    distsum, num = _distsum(sys.hamiltonian.intra, sys)
    for inter in sys.hamiltonian.inters
        d = _distsum(inter, sys)
        distsum += d[1]
        num += d[2]
    end
    return distsum/num
end

function _distsum(block, sys)
    distsum = 0.0
    num = 0
    for ((s1, s2), (target, source), _) in BlockIterator(block, sys.sysinfo)
        distsum += norm(site(sys, s2, source) - (site(sys, s1, target) + bravaismatrix(sys) * block.ndist))
        num += 1
    end
    return distsum, num
end

function uniquelinks(block::Block{Tv}, sys::System{E,L,T}) where {Tv,T,E,L}
    rdrs = [Tuple{SVector{E,T}, SVector{E,T}}[] for i in 1:nsublats(sys), j in 1:nsublats(sys)]
    bravais = bravaismatrix(sys)
    for (i, ((s1, s2), (target, source), (row, col), _)) in enumerate(BlockIterator(block, sys.sysinfo))
        if row > col || !iszero(block.ndist)
            rdr = _rdr(site(sys, s2, source), site(sys, s1, target) + bravais * block.ndist)
            push!(rdrs[s1, s2], rdr)
        end
    end
    return rdrs
end