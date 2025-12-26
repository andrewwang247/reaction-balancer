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


def display_solution(coefs: ndarray, mols: Tuple[str, ...]) -> str:
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
    left_mols = [parse(mol) for mol in left]
    right_mols = [parse(mol) for mol in right]
    if verbose:
        print('Molecules (L):')
        for mol, elem_counts in zip(left, left_mols):
            print(f'\t{mol}:', dict(elem_counts))
        print('Molecules (R):')
        for mol, elem_counts in zip(right, right_mols):
            print(f'\t{mol}:', dict(elem_counts))
    solutions = list(balance(left_mols, right_mols, verbose))
    if len(solutions) == 0:
        print('No solutions found.')
    else:
        print(f'{len(solutions)} solution{"" if len(solutions) == 1 else "s"}:')
    for left_coef, right_coef in solutions:
        left_disp = display_solution(left_coef, left)
        right_disp = display_solution(right_coef, right)
        print(f'\t{left_disp} -> {right_disp}')


if __name__ == '__main__':
    main()
