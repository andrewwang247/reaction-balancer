"""
Chemical reaction balance solver.

Copyright 2026. Andrew Wang.
"""
import logging
from typing import DefaultDict, Iterable, List, Tuple
from itertools import chain
import numpy as np
from sympy import Matrix, Rational  # type: ignore
from parse import parse

logger = logging.getLogger(__name__)


def __distinct_elems(mols: List[DefaultDict[str, int]]) -> List[str]:
    """Get the distinct elements that form the molecules."""
    elems = set()
    for mol in mols:
        for key in mol:
            elems.add(key)
    logger.info('Distinct elements (%d): %s', len(elems), elems)
    return list(elems)


def __scale_to_integers(rationals: List[Rational]) -> np.ndarray:
    """Scale a list of rationals to integers."""
    rational_rep = [num.as_numer_denom() for num in rationals]
    numers = np.array([rat[0] for rat in rational_rep])
    denoms = np.array([rat[1] for rat in rational_rep])
    lcm = np.lcm.reduce(denoms)
    coefs = numers * lcm / denoms
    coefs /= np.gcd.reduce(coefs)
    return coefs


def balance(left: List[DefaultDict[str, int]],
            right: List[DefaultDict[str, int]]) \
        -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """Balance the parsed left and ride sides."""
    elems = __distinct_elems(left + right)

    lin_sys = np.zeros((len(elems), len(left) + len(right)), dtype=int)
    for idx_elem, elem in enumerate(elems):
        for idx_mol, mol in enumerate(chain(left, right)):
            lin_sys[idx_elem, idx_mol] = mol[elem]
    lin_sys[:, len(left):] *= -1

    logger.info('Linear system of equations matrix:\n%s', lin_sys)
    nullspace: List[Matrix] = Matrix(lin_sys).nullspace(simplify=True)
    assert all(null_basis.shape[1] == 1 for null_basis in nullspace), \
        'Kernel basis should consist of column vectors.'
    kernel: List[List[Rational]] = [null_basis.flat()
                                    for null_basis in nullspace]

    logger.info('Nullity = %d', len(kernel))
    for ker in kernel:
        coefs = __scale_to_integers(ker)
        logger.info('Kernel basis vector %s scaled to %s', ker, coefs)
        if np.any(coefs < 0) and np.any(coefs > 0):
            continue
        coefs = np.abs(coefs)
        yield coefs[:len(left)], coefs[len(left):]


def solve(left: Iterable[str], right: Iterable[str]
          ) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """Balance left and right sides of chemical equation."""
    left_mols = [parse(mol) for mol in left]
    right_mols = [parse(mol) for mol in right]
    logger.info('Molecules (L): %s', list(left))
    logger.info('Molecules (R): %s', list(right))
    return balance(left_mols, right_mols)
