# Compute dimensions of cyclotomic (weighted) KLR(W)

This is [SageMath](https://www.sagemath.org/SageMath) code to compute dimensions of cyclotomic KLRW algebras. The code is due to Andrew Mathas and the original file can be found here: [Click](https://github.dev/AndrewAtLarge/GradedDimKLR). The purpose of this side is to give another explanation of how the code works and is also a (very mildly) updated version of the code.

# Background

The starting point is an amating paper by Hu and Shi [Click](https://arxiv.org/abs/2108.05508G) that proves the following dimension formula for cyclotomic KLR algebras:

```math
    \dim_q e(i)R^\Lambda_n e(j)
          = \sum_{w\in S_{i,j}}
              \prod_{t=1}^n[N^\Lambda(w,i,t)]_{i_t} q_{i_t}^{N^\Lambda(i,t)-1}
```
Here $R^\Lambda_n$ denotes the cyclotomic KLR algebra with $n$-strands, and $e(i)$ and $e(j)$ are idempotents corresponding to residue sequences. The sum runs over the subset of the symmteric group $S_n$ of permuations from $i$ to $j$. The $q$ indicates the grading and the $[N^\Lambda(w,i,t)]_{i_t}$ are certain quantum numbers that can be negative (careful with cancellations).

# What the code does

The code that can be downloaded here (klrw-dim.pv) provides a single command ``klr_cyclotomic_dimension`` that uses the Hu-Shi formula to compute the dimension. The syntax of this command is:

```python
   klr_cyclotomic_dimension(C, L, bi, bj=None, base=[], verbose=False)
```

where:

* ``C`` specifies the *Cartan type* of the quiver. This can either be a list
  of the form ``['A', 3]`` (finite type $A_3$), ``['A', 3, 1]`` (affine type
  $A_3^{(1)}$), ``['B',4]`` (finite type $B_4$), etc or ``C`` can be a [SageMath](https://www.sagemath.org/SageMath) Cartan type, such as ``CartanType(['D',3,1])``.

* ``L`` is a *list* that specifies the dominant weight. For example,
  ``[0,2,2,3]`` represents the dominant weight
  $\Lambda_0+2\Lambda_2+\Lambda_3$.

* ``i`` is an `n`-tuple of vertices of the quiver

* ``j`` is an `n`-tuple of vertices of the quiver. If ``j`` is omitted then
  ``j`` is set equal to ``i``.

* ``verbose`` an optional parameter that when set to ``True`` causes additional
  information from the Hu-Shi formula to be printed.

