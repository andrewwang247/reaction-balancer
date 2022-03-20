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
$ python main.py -l KNO3 -l C12H22O11 -r N2 -r CO2 -r H2O -r K2CO3 -v
Left side:  KNO3, C12H22O11
Right side:  N2, CO2, H2O, K2CO3
Elements left: {'K': 1, 'N': 1, 'O': 3}, {'C': 12, 'H': 22, 'O': 11}
Elements right: {'N': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'K': 2, 'C': 1, 'O': 3}
Distinct elements: N, C, O, K, H
Linear system of equations matrix
[[ 1  0 -2  0  0  0]
 [ 0 12  0 -1  0 -1]
 [ 3 11  0 -2 -1 -3]
 [ 1  0  0  0  0 -2]
 [ 0 22  0  0 -2  0]]
Nullity = 1
Kernel basis Matrix([[2], [5/24], [1], [3/2], [55/24], [1]]) simplified to [48 5 24 36 55 24]
48 KNO3 + 5 C12H22O11 -> 24 N2 + 36 CO2 + 55 H2O + 24 K2CO3

$ python3 main.py -l "Al(OH)3" -l H2SO4 -r "Al2(SO4)3" -r H2O -v
Left side:  Al(OH)3, H2SO4
Right side:  Al2(SO4)3, H2O
Elements left: {'Al': 1, 'O': 3, 'H': 3}, {'H': 2, 'S': 1, 'O': 4}
Elements right: {'Al': 2, 'S': 3, 'O': 12}, {'H': 2, 'O': 1}
Distinct elements: O, Al, H, S
Linear system of equations matrix
[[  3   4 -12  -1]
 [  1   0  -2   0]
 [  3   2   0  -2]
 [  0   1  -3   0]]
Nullity = 1
Kernel basis Matrix([[1/3], [1/2], [1/6], [1]]) simplified to [2 3 1 6]
2 Al(OH)3 + 3 H2SO4 -> Al2(SO4)3 + 6 H2O
```

Execute unit tests for parsing with `pytest`.

## Methodology

For each chemical reaction, the balancer counts the elements in each component. It uses these values to generate a system of linear equations that represent the constraints of a possible solution. A basis of the corresponding matrix null space provides the coefficients of the balanced equation.

An example is shown below:

```text
Problem: a CH4 + b O2 -> c CO2 + d H2O

C constraint: 1a + 0b = 1c + 0d
O constraint: 0a + 2b = 2c + 1d
H constraint: 4a + 0b = 0c + 2d

Matrix form:
[[ 1  0 -1  0]
 [ 0  2 -2 -1]
 [ 4  0  0 -2]]

Kernel vector: [1 2 1 2]

Solution: CH4 + 2 O2 -> CO2 + 2 H2O
```
