"""
Balance user-provided chemical reactions.

Copyright 2026. Andrew Wang.
"""
import logging
from typing import Tuple
from click import command, option
import numpy as np
import numpy.typing as npt
from src import solve
# pylint: disable=no-value-for-parameter


def _display_solution(coefs: npt.NDArray[np.int_],
                      mols: Tuple[str, ...]) -> str:
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
    logging.basicConfig(level=logging.INFO if verbose else logging.WARN)
    solutions = list(solve(left, right))
    if len(solutions) == 0:
        print('No solutions found.')
    else:
        print(f'Solutions ({len(solutions)}):')
    for left_coef, right_coef in solutions:
        left_disp = _display_solution(left_coef, left)
        right_disp = _display_solution(right_coef, right)
        print(f'\t{left_disp} -> {right_disp}')


if __name__ == '__main__':
    main()
