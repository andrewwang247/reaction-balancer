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

The custom parser supports standard chemical formulas consisting of elements, subscripts, and nested parentheses. 

Some examples:

```text
$ python3 main.py -l KNO3 -l C12H22O11 -r N2 -r CO2 -r H2O -r K2CO3 -v
Molecules (L):
	KNO3: {'K': 1, 'N': 1, 'O': 3}
	C12H22O11: {'C': 12, 'H': 22, 'O': 11}
Molecules (R):
	N2: {'N': 2}
	CO2: {'C': 1, 'O': 2}
	H2O: {'H': 2, 'O': 1}
	K2CO3: {'K': 2, 'C': 1, 'O': 3}
Distinct elements (5): ['O', 'C', 'H', 'K', 'N']
Linear system of equations matrix:
[[ 3 11  0 -2 -1 -3]
 [ 0 12  0 -1  0 -1]
 [ 0 22  0  0 -2  0]
 [ 1  0  0  0  0 -2]
 [ 1  0 -2  0  0  0]]
Nullity = 1. Kernel basis:
	[2, 5/24, 1, 3/2, 55/24, 1]
Solutions (1):
	48 KNO3 + 5 C12H22O11 -> 24 N2 + 36 CO2 + 55 H2O + 24 K2CO3
```

```text
$ python3 main.py -l "Al(OH)3" -l H2SO4 -r "Al2(SO4)3" -r H2O -v
Molecules (L):
	Al(OH)3: {'Al': 1, 'O': 3, 'H': 3}
	H2SO4: {'H': 2, 'S': 1, 'O': 4}
Molecules (R):
	Al2(SO4)3: {'Al': 2, 'S': 3, 'O': 12}
	H2O: {'H': 2, 'O': 1}
Distinct elements (4): ['Al', 'O', 'H', 'S']
Linear system of equations matrix:
[[  1   0  -2   0]
 [  3   4 -12  -1]
 [  3   2   0  -2]
 [  0   1  -3   0]]
Nullity = 1. Kernel basis:
	[1/3, 1/2, 1/6, 1]
Solutions (1):
	2 Al(OH)3 + 3 H2SO4 -> Al2(SO4)3 + 6 H2O
```

```text
$ python3 main.py -l H -l P -l O -r H -r P -r O2 -v
Molecules (L):
	H: {'H': 1}
	P: {'P': 1}
	O: {'O': 1}
Molecules (R):
	H: {'H': 1}
	P: {'P': 1}
	O2: {'O': 2}
Distinct elements (3): ['O', 'H', 'P']
Linear system of equations matrix:
[[ 0  0  1  0  0 -2]
 [ 1  0  0 -1  0  0]
 [ 0  1  0  0 -1  0]]
Nullity = 3. Kernel basis:
	[1, 0, 0, 1, 0, 0]
	[0, 1, 0, 0, 1, 0]
	[0, 0, 2, 0, 0, 1]
Solutions (3):
	H + 0 P + 0 O -> H + 0 P + 0 O2
	0 H + P + 0 O -> 0 H + P + 0 O2
	0 H + 0 P + 2 O -> 0 H + 0 P + O2
```

```text
$ python3 main.py -l C -r Ne -v
Molecules (L):
	C: {'C': 1}
Molecules (R):
	Ne: {'Ne': 1}
Distinct elements (2): ['Ne', 'C']
Linear system of equations matrix:
[[ 0 -1]
 [ 1  0]]
Nullity = 0. Kernel basis:
No solutions found.
```

## Methodology

For each chemical reaction, the balancer counts the elements in each component. It uses these values to generate a system of linear equations that represent the constraints of a possible solution. Turning the system into a matrix and computing the null space yields solutions for the coefficients of the balanced equation.

We work through an example for $a \cdot CH_4 + b \cdot O_2 \to c \cdot CO_2 + d \cdot H_2 O$.
- $C$ constraint: $1a + 0b = 1c + 0d$
- $O$ constraint: $0a + 2b = 2c + 1d$
- $H$ constraint: $4a + 0b = 0c + 2d$

Computing the matrix null space
```math
\ker
\begin{pmatrix}
1 & 0 & -1 & 0 \\
0 & 2 & -2 & -1 \\
4 & 0 & 0 & -2
\end{pmatrix}
= \text{span}
\begin{pmatrix}
1 & 2 & 1 & 2
\end{pmatrix}
```

Thus, the balanced equation is $CH_4 + 2 \\, O_2 \to CO_2 + 2 \\, H_2 O$.

## Testing

Both parsing and balancing are validated with extensive unit testing. The test cases for parsing and balancing are stored in `tst/molecules.json` and `tst/equations.json`, respectively. Run `pytest` to run the full suite of tests.
