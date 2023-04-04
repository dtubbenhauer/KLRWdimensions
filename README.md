# Compute dimensions of cyclotomic (weighted) KLR(W)

This is [SageMath](https://www.sagemath.org/SageMath) code to compute dimensions of cyclotomic KLRW algebras. The code is due to Andrew Mathas (**all the credits should go to Andrew!**) and the original file and explanation can be found here: [click](https://github.dev/AndrewAtLarge/GradedDimKLR). The purpose of this side is to give another explanation of how the code works and is also a (very mildly) updated version of the code.

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

Before we do examples, let me explain how to get started.

# Getting started

* Using the online calculator [SageMathCell](https://sagecell.sagemath.org/) is easy: Cut-and-paste the text from *klrw-dim.py* into a cell
  and then type the ``klr_cyclotomic_dimension`` commands into the bottom of the cell.

```python
   klr_cyclotomic_dimension(['A',3],[2,3], [2,3,3,2,1], [2,3,2,3,1])
```
![Screenshot of the online calculator and the code.](https://github.com/dtubbenhauer/KLRWdimensions/blob/main/sagemath-online-calc.png)

* With a local installation of [SageMath](https://www.sagemath.org/), start [SageMath](https://www.sagemath.org/) and attach the file ``graded_dim_klr`` using:

```python
   sage: %attach graded_dim_klr
   sage: klr_cyclotomic_dimension(['A',3],[2,3], [2,3,3,2,1], [2,3,2,3,1])
```
# Examples

Ok, let us have a look at
```python
   klr_cyclotomic_dimension(['A',3],[2,3], [2,3,3,2,1], [2,3,2,3,1])
```
The output is
```python
   ({(q + 1/q)^2*q: [(0,2,1,3)]}, (q^2 + 1)^2/q)
```
This is saying that there is only one diagram with bottom sequence 2,3,3,2,1 and top sequence 2,3,2,3,1, namely:

![KLRW example.](https://github.com/dtubbenhauer/KLRWdimensions/blob/main/klrw-diagram.png)

Note that I read from right to left. The bottom sequence of numbers is the position of the strings. The degree of the diagram is (q^2 + 1)^2/q.permuation

All other outputs are read similarly. Note that we can get more than one diagram, for example
```python
   klr_cyclotomic_dimension(['B',3],[2,3], [2,3,3,2,1], [2,3,2,3,1])
```
gives
```python
   ({(q^2 + 1/q^2)*q^4: [(0,2,3)],
  (q^2 + 1/q^2 + 1)*(q^2 + 1/q^2)*q^4: [(0,2,1,3)]},
 (q^4 + 1)*(q^2 + 1)^2)
```
Thus, there are two relevant diagrams, one determined by the permutation (0,2,3) and the other by (0,2,1,3).

# What the code can do - more advanced

If we turn on ``verbose=True``, then for example
```python
   klr_cyclotomic_dimension(['B',3],[2,3], [2,3,3,2,1], [2,3,2,3,1], verbose=True)
```
we get
```python
Subgroup of permutations = [(2,3), (0,2,3), (1,3,2), (0,2,1,3)]
N(1,t)-1: 0  2  0  0  1
N(w,t):   1  3  3  0  2
N(w,t):   1  1  1  1  2
N((0,2,3),t):   1  1  1  1  2
X((0,2,3)): (q^4 + 1)*q^2
N(w,t):   1  3  1  0  2
N(w,t):   1  3  1  1  2
N((0,2,1,3),t):   1  3  1  1  2
X((0,2,1,3)): (q^4 + 1)*(q^2 + q + 1)*(q^2 - q + 1)
(q^2 + 1/q^2)*q^4: (0,2,3)
(q^2 + 1/q^2 + 1)*(q^2 + 1/q^2)*q^4: (0,2,1,3)

({(q^2 + 1/q^2)*q^4: [(0,2,3)],
  (q^2 + 1/q^2 + 1)*(q^2 + 1/q^2)*q^4: [(0,2,1,3)]},
 (q^4 + 1)*(q^2 + 1)^2)
```
which is saying that the code tried four permutations, but only two are relevant. In this example there are no cancellations, but try for example
```python
klr_cyclotomic_dimension(['F',4],[2,3], [2,3,3,2,1,2,3,3,2], [2,3,2,3,1,2,3,3,2], verbose=True)
```
to see cancellations.

Turning on the base, for example
```python
   klr_cyclotomic_dimension(['B',3],[2,3], [3,3,2,1], [3,2,3,1], base=[2])
```
makes the computation faster, but might get the wrong result. In the above example we get
```python
   ({}, 0)
```
which is wrong. We should get
```python
   ({(q^2 + 1/q^2)*q^4: [(0,2,3)],
  (q^2 + 1/q^2 + 1)*(q^2 + 1/q^2)*q^4: [(0,2,1,3)]},
 (q^4 + 1)*(q^2 + 1)^2)
```
but since we excluded the first entry, the 0th position, the above permutations cannot appear, so we get zero as the result.

``base`` is still useful for larger computations, but **use with care**.

# Contact

If there are any questions, then please feel free to contact me: dtubbenhauer@gmail.com

Have fun playing with the code!
