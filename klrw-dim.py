r"""
    Compute the graded dimension of  e(i)R^Lambda_\alpha e(j),
    where \Lambda is a dominant weight and i,j\in I^\alpha.

    By Hu-Shi, if i,j\in I^n then

    .. math:

        \dim_q e(i)R^\Lambda_\alpha e(j)
          = \sum_{w\in S_{(i,j)}}
              \prod_{t=1}^n[N^\Lambda(w,i,t)]_{i_t} q_{i_i}^{N^\Lambda(i,t)-1}}

    sage: R = RootSystem(['G', 2])
    sage:  klr_cyclotomic_dimension(R, [1], [1,2,1,1,2])

    sage: KLRW_Idempotents(['F',4],[2],6,latex=True)


    Andrew Mathas
"""

from sage.combinat.combination             import Combinations
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.tuple                   import Tuples
from sage.graphs.graph_plot                import GraphPlot
from sage.groups.perm_gps.permgroup_named  import SymmetricGroup
from sage.plot.colors                      import Color
from sage.rings.infinity                   import Infinity
from sage.rings.real_mpfr                  import RR
from sage.structure.sage_object            import SageObject
from sage.typeset.ascii_art                import AsciiArt

from decimal import Decimal

def quantum_integer(q,k):
    '''
    Return the quantum integer (q^k-q^{-k})/(q-q^{-1})

        [1] = 1
        [2] = q + q^-1
        [3] = q^2 + 1 + q-2
        [-k] = -[k]

    '''
    if k < 0:
        return -quantum_integer(q, -k)

    return sum( q**(k-2*i-1) for i in range(k))


def klr_cyclotomic_dimension(C, L, bi, bj=None, base=[], verbose=False):
    r"""
    This function implements the Hu-Shi formula for the graded dimension of the
    weight space e(i) R^L_n e(j), where i,j \in I^n. Here:

        - `C`  is the Cartan type
        - `L`  is a list specifying the dominant weight
        - `bi` is a weight in I^n
        - `bj` is a weight in I^n, which defaults to `bi`

    If `verbose` is set to `True` then some details of the computation are
    printed as well.

    EXAMPLES::

        sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1])
        Subgroup of permutations = < () >
        N(1,t)-1: 0  0  0  0
        N(w,t):   1  1  1  1
        X(w): 1
        1
        sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1], [2,3,4,1])
        Subgroup of permutations = < () >
        N(1,t)-1: 0  0  0  0
        N(w,t):   1  1  1  1
        X(w): 1
        1
        sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1], [2,4,3,1])
        Subgroup of permutations = < (1,2) >
        N(1,t)-1: 0  0  0  0
        N(w,t):   1  1  1  1
        X(w): 1
        1
        sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,1,2])
        Subgroup of permutations = < () >
        N(1,t)-1: 0  0  1
        N(w,t):   1  1  2
        X(w): q^2 + 1
        q^2 + 1
        sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,2,1])
        Subgroup of permutations = < () >
        N(1,t)-1: 0  0  1
        N(w,t):   1  1  2
        X(w): q^2 + 1
        q^2 + 1
        sage:
        sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,1,2], [0,2,1])
        Subgroup of permutations = < (1,2) >
        N(1,t)-1: 0  0  1
        N(w,t):   1  1  1
        X(w): q
        q
        sage: klr_cyclotomic_dimension(['A',1],[1,1],[1],[1])
        Subgroup of permutations = < () >
        N(1,t)-1: 1
        N(w,t):   2
        X(w): q^2 + 1
        q^2 + 1
        sage: klr_cyclotomic_dimension(['A',1],[1,1],[1],[1])
        Subgroup of permutations = < () >
        N(1,t)-1: 1
        N(w,t):   2
        X(w): q^2 + 1
        q^2 + 1
        sage: klr_cyclotomic_dimension(['B',3],[2], [2,3,3,2,1])
        Subgroup of permutations = < (), (1,2), (0,3), (0,3)(1,2) >
        N(1,t)-1: 0  1 -1  0  1
        N(w,t):   1  2  0  1  2
        X(w): 0
        N(w,t):   1  2  2  1  2
        X(w): (q^4 + 1)*(q^2 + 1)^2/q^2
        N(w,t):   1  0 -2  1  2
        X(w): 0
        N(w,t):   1  0  0  1  2
        X(w): 0
        (q^4 + 1)*(q^2 + 1)^2/q^2
    """
    if verbose:
        def vprint(*args): print(' '.join(f'{a}' for a in args))
    else:
        def vprint(*args): pass

    # bj defaults to bi
    if bj is None:
        bj = bi

    if len(bi) != len(bj):
        raise TypeError(f'{bi} and {bj} must have the same length')

    # Convert C to Cartantype
    try:
        C = CartanType(C)
    except TypeError:
        raise TypeError(f'{C} must be a Cartan type')

    # index set for quiver
    I = C.index_set()
    #vprint(f'{I=}')

    # if bi and bj are strings convert them to a list
    if isinstance(bi,str):
        bi = [Integer(i) for i in bi]

    if isinstance(bj,str):
        bj = [Integer(j) for j in bj]

    if isinstance(base,str):
        base = [Integer(j) for j in base]

    if any(i not in I for i in bi):
        raise TypeError(f'{bi} must in I^n for {I=}')

    if any(j not in I for j in bj):
        raise TypeError(f'{bj} must in I^n for {I=}')

    # shorthand for Cartan matrix entries
    Cij = lambda i,j: C.cartan_matrix()[I.index(i), I.index(j)]

    # shorthand for q_di, where is the Cartan symmetriser
    qdi  = lambda i: q**(C.symmetrizer()[i])

    # work in SymmetricGroup(n)
    n = len(bi)

    q = var('q')

    N = Integer(0) # will become \sum_w \prod_t [N^L(w,bi,t)]q^{N^L(1,bi,t)-1}

    Sn = SymmetricGroup(range(n))
    SN = SymmetricGroup(range(n+len(base)))

    # find the subgroup of Sn that maps bi to bj
    generators = []
    for w in Sn:
        if all( bj[w(i)]==bi[i] for i in range(n) ) and not w.is_one():
            generators.append( prod(SN.simple_reflection(i+len(base)) for i in w.reduced_word()) )

    i = 0
    while i < len(base)-1:
        j = i+1
        while j < len(base) and Cij(base[i], base[j])==0:
            j += 1
        if j<len(base) and base[i] == base[j]:
            generators.append( SN((i,j)) )
        i += 1

    bi = base + bi
    bj = base + bj
    n = len(bi)
    Nwt = lambda w,t: L.count(bi[t]) - sum(Cij(bi[t], bi[j]) for j in range(t) if w(j)<w(t))

    vprint(f'Subgroup of permutations = {generators}') #.replace('[', '< ').replace(']',' >'))
    # first calculate NOne[t] = Nwt(1, t) for t in range(n)
    tally = {}
    NOne = [ Nwt(SN.one(), t)-1 for t in range(n)]
    vprint(f'N(1,t)-1:{" ".join(f"{i:>2}" for i in NOne)}')
    #vprint('qdi:   ',' '.join(f'{qdi(bi[t])}' for t in range(n)))
    subgroup = SN.subgroup(generators)
    for w in subgroup:
        # Sym(i,j) = { w\in Sym_n | wi = j}
        if all( bj[w(i)]==bi[i] for i in range(n) ): # if w(bi) == bj
            vprint(f'N(w,t): ',' '.join(f'{Nwt(w,t):>2}' for t in range(n)))
            Nw = prod([quantum_integer( qdi(bi[t]), Nwt(w,t) )*( qdi(bi[t])**(NOne[t]) ) for t in range(n) ])
            if Nw != 0:
                vprint(f'N({w},t): ',' '.join(f'{Nwt(w,t):>2}' for t in range(n)))
                vprint(f'X({w}):', Nw.factor())
                if -Nw in tally:
                    tally[-Nw].pop()
                    if tally[-Nw] == []:
                        # tally is sorted by length so we remove a w of longest length
                        del tally[-Nw] 
                else:
                    if Nw not in tally:
                        tally[Nw] = []
                    tally[Nw].append(SN(w))
                    tally[Nw].sort(key=lambda w: w.length())
                N += Nw

    for Nw in tally:
        vprint(f'{Nw}: {", ".join(f"{w}" for w in tally[Nw])}')

    if N != 0:
        N = N.factor()

    return tally, N
