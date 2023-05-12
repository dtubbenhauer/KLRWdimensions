r"""
The following code runs in SageMath, e.g. in the online version available under https://sagecell.sagemath.org/

The main definition is klr_cyclotomic_dimension(C, L, bi, bj, verbose=True). This can be used as follows:

- the first entry is the type, e.g. ['A',4]. The code below actually computes the general formula of Hu-Shi and the type can be any input, but for the present work type A is needed and the second entry is the number of different F_i's used for the F forms

- L is the level in terms of fundamental weights. For example, to get (4,4,4,0,0,0,0) use [3,3,3,3] so the 3rd fundamental weight (1,1,1,0,0,0,0) four times

- bi and bj are the exploded residue sequences of the webs in question

- the final entry set to false suppresses all intermediate outputs and only gives the main result

- the code computes KLR dimensions, so e.g. the exploded and not the non-exploded residue sequences

"""

r"""
Some preliminaries loaded
"""

from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.tuple                   import Tuples
from sage.graphs.graph_plot                import GraphPlot
from sage.groups.perm_gps.permgroup_named  import SymmetricGroup
from sage.plot.colors                      import Color
from sage.rings.real_mpfr                  import RR
from sage.structure.sage_object            import SageObject
from sage.typeset.ascii_art                import AsciiArt

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
    
def klr_cyclotomic_dimension(C, L, bi, bj, verbose=True):
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
        1
        sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1], [2,4,3,1])
        1
        sage: sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,1,2])
        q^2 + 1
        sage: sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,2,1])
        q^2 + 1
    sage: sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,1,2], [0,2,1])
    q

    """
    if verbose:
        def vprint(*args): print(' '.join(f'{a}' for a in args))
    else:
        def vprint(*args): pass

    # bj defaults to bi
    bj = bj

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

    N = 0 # will become \sum_w \prod_t [N^L(w,bi,t)]q^{N^L(1,bi,t)-1}

    Nwt = lambda w,t: L.count(bi[t]) - sum(Cij(bi[t], bi[j]) for j in range(t) if w(j)<w(t))
    Sn = SymmetricGroup(range(n))
    # find the subset of Sn that maps bi to bj
    generators = []
    positions={}
    for w in Sn:
        if all( bj[w(i)]==bi[i] for i in range(n) ):
            generators.append(w)
    vprint(f'Subset of permutations=', generators)
    # first calculate NOne[t] = Nwt(1, t) for t in range(n) 
    NOne = [ Nwt(Sn.one(), t)-1 for t in range(n)]
    vprint(f'N(unit,t)-1:  {" ".join(f"{i}" for i in NOne)}')
    #vprint('qdi:   ',' '.join(f'{qdi(bi[t])}' for t in range(n)))
    for w in Sn.subgroup( generators ):
        # Sym(i,j) = { w\in Sym_n | wi = j}
        if all( bj[w(i)]==bi[i] for i in range(n) ): # if w(bi) == bj
            vprint(f'N(w,t):  {w} '+' '.join(f'{Nwt(w,t)}' for t in range(n)))
            Nw = prod([
                quantum_integer( qdi(bi[t]), Nwt(w,t) )*( qdi(bi[t])**(NOne[t]) )
                for t in range(n)
            ])
            vprint(f'X(w):', Nw.expand())
            N += Nw
    return N.expand()

r"""

The following code can be run to verify the examples in the note, but can also be adjusted to any example the reader wants to run.

The first two are Example 4.5(a), the second is Example 4.5(b). Remove the # to run

The final four are the four summands in Example 5.7. The used symmetric group in this case, so the calculation might take a while

"""
#klr_cyclotomic_dimension(['A',1],[1,1],[1],[1])
#klr_cyclotomic_dimension(['A',1],[1,1],[1,1],[1,1])
#klr_cyclotomic_dimension(['A',3],[2,2,2],[2,3,2,2,1],[2,3,1,2,2])

#klr_cyclotomic_dimension(['A',5],[2,2],[4,5,3,4,2,3,1,2,4,3,5,4,2,3,1,2],[2,2,3,3,4,4,5,5,1,1,2,2,3,3,4,4])
#klr_cyclotomic_dimension(['A',5],[2,2],[4,5,3,4,2,3,2,1,4,3,5,4,2,3,1,2],[2,2,3,3,4,4,5,5,1,1,2,2,3,3,4,4])
#klr_cyclotomic_dimension(['A',5],[2,2],[4,5,3,4,3,2,1,2,4,3,5,4,2,3,1,2],[2,2,3,3,4,4,5,5,1,1,2,2,3,3,4,4])
#klr_cyclotomic_dimension(['A',5],[2,2],[4,5,3,4,3,2,2,1,4,3,5,4,2,3,1,2],[2,2,3,3,4,4,5,5,1,1,2,2,3,3,4,4])
