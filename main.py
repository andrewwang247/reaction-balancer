"""
Balance user-provided chemical reactions.

Copyright 2022. Andrew Wang.
"""
from typing import Tuple
from click import command, option
from numpy import ndarray
from parse import parse
from balance import balance
# pylint: disable=no-value-for-parameter


def to_string(coefs: ndarray, mols: Tuple[str, ...]) -> str:
    """Construct string of coefficients and molecules."""
    multipled = []
    for coef, mol in zip(coefs, mols):
        multipled.append(f'{coef} {mol}' if coef != 1 else mol)
    return ' + '.join(multipled)


@command()
@option('--left', '-l', type=str, required=True, multiple=True,
        help='Left side of chemical reaction.')
@option('--right', '-r', type=str, required=True, multiple=True,
        help='Right side of chemical reaction.')
@option('--verbose', '-v', is_flag=True, default=False,
        help='Set verbosity of solving process.')
def main(left: Tuple[str, ...], right: Tuple[str, ...], verbose: bool):
    """Balance user-provided chemical reactions."""
    if verbose:
        print('Molecules (L): ', ', '.join(left))
        print('Molecules (R): ', ', '.join(right))
    left_mols = [parse(mol) for mol in left]
    right_mols = [parse(mol) for mol in right]
    if verbose:
        print('Elements (L):', ', '.join(
            map(str, (dict(mol) for mol in left_mols))))
        print('Elements (R):', ', '.join(
            map(str, (dict(mol) for mol in left_mols))))
    solutions = balance(left_mols, right_mols, verbose)
    for left_coef, right_coef in solutions:
        left_disp = to_string(left_coef, left)
        right_disp = to_string(right_coef, right)
        print(left_disp, '->', right_disp)


if __name__ == '__main__':
    main()
