# Reaction Balancer

A chemical reaction balancer that solves user-provided equations.

```text
Usage: main.py [OPTIONS]

  Balance user-provided chemical reactions.

Options:
  -l, --left TEXT   Left side of chemical reaction.  [required]
  -r, --right TEXT  Right side of chemical reaction.  [required]
  -v, --verbose     Set verbosity of solving process.
  --help            Show this message and exit.
```

## Grammar

The custom parser supports standard chemical formulas consisting of elements, subscripts, and nested parentheses. Some examples:

```text
$ python3 main.py -l KNO3 -l C12H22O11 -r N2 -r CO2 -r H2O -r K2CO3 -v
INFO:balance:Molecules (L): ['KNO3', 'C12H22O11']
INFO:balance:Molecules (R): ['N2', 'CO2', 'H2O', 'K2CO3']
INFO:balance:Distinct elements (5): {'N', 'H', 'O', 'K', 'C'}
INFO:balance:Linear system of equations matrix:
[[ 1  0 -2  0  0  0]
 [ 0 22  0  0 -2  0]
 [ 3 11  0 -2 -1 -3]
 [ 1  0  0  0  0 -2]
 [ 0 12  0 -1  0 -1]]
INFO:balance:Nullity = 1
INFO:balance:Kernel basis vector [2, 5/24, 1, 3/2, 55/24, 1] scaled to [48 5 24 36 55 24]
Solutions (1):
    48 KNO3 + 5 C12H22O11 -> 24 N2 + 36 CO2 + 55 H2O + 24 K2CO3
```

```text
$ python3 main.py -l "Al(OH)3" -l H2SO4 -r "Al2(SO4)3" -r H2O -v
INFO:balance:Molecules (L): ['Al(OH)3', 'H2SO4']
INFO:balance:Molecules (R): ['Al2(SO4)3', 'H2O']
INFO:balance:Distinct elements (4): {'S', 'H', 'Al', 'O'}
INFO:balance:Linear system of equations matrix:
[[  0   1  -3   0]
 [  3   2   0  -2]
 [  1   0  -2   0]
 [  3   4 -12  -1]]
INFO:balance:Nullity = 1
INFO:balance:Kernel basis vector [1/3, 1/2, 1/6, 1] scaled to [2 3 1 6]
Solutions (1):
    2 Al(OH)3 + 3 H2SO4 -> Al2(SO4)3 + 6 H2O
```

```text
$ python3 main.py -l H -l P -l O -r H -r P -r O2 -v
INFO:balance:Molecules (L): ['H', 'P', 'O']
INFO:balance:Molecules (R): ['H', 'P', 'O2']
INFO:balance:Distinct elements (3): {'P', 'O', 'H'}
INFO:balance:Linear system of equations matrix:
[[ 0  1  0  0 -1  0]
 [ 0  0  1  0  0 -2]
 [ 1  0  0 -1  0  0]]
INFO:balance:Nullity = 3
INFO:balance:Kernel basis vector [1, 0, 0, 1, 0, 0] scaled to [1 0 0 1 0 0]
INFO:balance:Kernel basis vector [0, 1, 0, 0, 1, 0] scaled to [0 1 0 0 1 0]
INFO:balance:Kernel basis vector [0, 0, 2, 0, 0, 1] scaled to [0 0 2 0 0 1]
Solutions (3):
    H + 0 P + 0 O -> H + 0 P + 0 O2
    0 H + P + 0 O -> 0 H + P + 0 O2
    0 H + 0 P + 2 O -> 0 H + 0 P + O2
```

```text
$ python3 main.py -l C -r Ne -v
INFO:balance:Molecules (L): ['C']
INFO:balance:Molecules (R): ['Ne']
INFO:balance:Distinct elements (2): {'Ne', 'C'}
INFO:balance:Linear system of equations matrix:
[[ 0 -1]
 [ 1  0]]
INFO:balance:Nullity = 0
No solutions found.
```

## Methodology

For each chemical reaction, the balancer counts the elements in each component. It uses these values to generate a system of linear equations that represent the constraints of a possible solution. Turning the system into a matrix and computing the null space yields solutions for the coefficients of the balanced equation. We work through an example for $a \cdot CH_4 + b \cdot O_2 \to c \cdot CO_2 + d \cdot H_2 O$. Each element constrains the coefficients in the form of a linear equation.

- $C$ constraint: $1a + 0b = 1c + 0d$ or $1a + 0b - 1c - 0d = 0$
- $O$ constraint: $0a + 2b = 2c + 1d$ or $0a + 2b - 2c - 1d = 0$
- $H$ constraint: $4a + 0b = 0c + 2d$ or $4a + 0b - 0c - 2d = 0$

In matrix form, this system of linear equations is:

```math
\begin{pmatrix}
1 & 0 & -1 & 0 \\
0 & 2 & -2 & -1 \\
4 & 0 & 0 & -2
\end{pmatrix}
\begin{pmatrix}
a \\ b \\ c \\ d
\end{pmatrix}
=
\vec{0}
```

We compute the null space and represent it as a set of basis vectors with rational components. These vectors are the unique solutions for the chemical coefficients. We rescale such that each component is an integer before presenting the solution.

```math
\ker
\begin{pmatrix}
1 & 0 & -1 & 0 \\
0 & 2 & -2 & -1 \\
4 & 0 & 0 & -2
\end{pmatrix}
= \text{span}
\begin{pmatrix}
1 \\ 2 \\ 1 \\ 2
\end{pmatrix}
```

Thus, the balanced equation is $CH_4 + 2 \\, O_2 \to CO_2 + 2 \\, H_2 O$.

## Testing

Both parsing and balancing are validated with extensive unit testing. The test cases for parsing and balancing are stored in `test/molecules.json` and `test/equations.json`, respectively. Run `pytest` to execute the full suite of tests.
