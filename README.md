# Compute dimensions of cyclotomic (weighted) KLR(W)

This is [SageMath](https://www.sagemath.org/SageMath) code to compute dimensions of cyclotomic KLRW algebras. The code is due to Andrew Mathas (**all the credits should go to Andrew!**) and the original file can be found here: [click](https://github.dev/AndrewAtLarge/GradedDimKLR). The purpose of this side is to give another explanation of how the code works and is also a (very mildly) updated version of the code.

# Background

The starting point is an amating paper by Hu and Shi [click](https://arxiv.org/abs/2108.05508G) that proves the following dimension formula for cyclotomic KLR algebras:

```math
    \dim_q e(j)R^\Lambda_n e(i)
          = \sum_{w\in S_{i,j}}
              \prod_{t=1}^n[N^\Lambda(w,i,t)]_{i_t} q_{i_t}^{N^\Lambda(i,t)-1}
```
Here $R^\Lambda_n$ denotes the cyclotomic KLR algebra with $n$-strands, and $e(i)$ and $e(j)$ are idempotents corresponding to residue sequences. The sum runs over the subset of the symmteric group $S_n$ of permuations from $i$ to $j$. The $q$ indicates the grading and the $[N^\Lambda(w,i,t)]_{i_t}$ are certain quantum numbers that can be negative (careful with cancellations).

# What the code can do - basics

The code that can be downloaded here (klrw-dim.pv) provides a single command ``klr_cyclotomic_dimension`` that uses the Hu-Shi formula to compute the dimension. The syntax of this command is:

```python
   klr_cyclotomic_dimension(C, L, bi, bj=None, base=[], verbose=False)
```

where:

* The first entry ``C`` specifies the **Cartan type** of the quiver that is used for the KLR algebra. The Crtan types are as in [CartanType]{https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/cartan_type.html} so for example ``['A', 3]`` (finite type $A_3$), ``['A', 3, 1]`` (affine type  $A_3^{(1)}$), ``['B',4]`` (finite type $B_4$), etc.

* ``L`` is a list that specifies the dominant weight, that is, the **number and the labels of the red strings** in the KLRW world. For example,
  ``[0,2,2,3]`` represents the dominant weight $\Lambda_0+2\Lambda_2+\Lambda_3$.

* ``i`` is an n-tuple of vertices of the quiver, the residue sequence = **coloring of the strings** at the, say, *bottom* of a diagram. (This does not need to be an n-tuple, see ``base`` below. Ditto for the next entry.)

* ``j`` is an n-tuple of vertices of the quiver, the residue sequence = **coloring of the strings** at the, say, *top* of a diagram. If ``j`` is omitted, then ``j`` is set equal to ``i``.

* ``base`` is a tuple that is similar to ``i`` and ``j`` and serves as the **base**. If this is not empty (it is empty by default), then the residue sequences are actually the concatenations of ``base`` and ``i`` respectively ``j``. The idea is that the Hu-Shi formula runs voer large symmetric groups, which might be outside of the calculation scope. With a nontrivial base the calculation of the permutations is restricted to ``i`` and ``j`` only which makes the calculation faster. This however **should be treated with care** as the Hu-Shi formual is not local.

* ``verbose`` an optional parameter that when set to ``True`` gives additional output that might be useful to analize what is going on.

