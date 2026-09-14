"""Chemical reaction balance solver.

Copyright 2026. Andrew Wang.
"""

import logging
from itertools import chain
from typing import TYPE_CHECKING

import numpy as np
from sympy import Matrix, Rational

if TYPE_CHECKING:
    from .parse import Elements

logger = logging.getLogger(__name__)

type IntArr = np.ndarray[tuple[int], np.dtype[np.int_]]
type ObjArr = np.ndarray[tuple[int], np.dtype[np.object_]]


def _distinct_elems(mols: list[Elements]) -> list[str]:
    """Get the distinct elements that form the molecules."""
    elems: set[str] = set()
    for mol in mols:
        for key in mol:
            elems.add(key)
    logger.info("Distinct elements (%d): %s", len(elems), elems)
    return list(elems)


def _scale_to_integers(rationals: list[Rational]) -> IntArr:
    """Scale a list of rationals to integers."""
    rational_rep = [num.as_numer_denom() for num in rationals]
    numers = np.array([rat[0] for rat in rational_rep])
    denoms = np.array([rat[1] for rat in rational_rep])
    lcm = np.lcm.reduce(denoms)
    coefs: ObjArr = numers * lcm / denoms
    coefs /= np.gcd.reduce(coefs)
    return coefs.astype(np.int_)


def solve(
    left: list[Elements],
    right: list[Elements],
) -> list[tuple[IntArr, IntArr]]:
    """Balance left and right sides of chemical equation."""
    elems = _distinct_elems(left + right)
    lin_sys = np.zeros((len(elems), len(left) + len(right)), dtype=np.int_)

    for idx_elem, elem in enumerate(elems):
        for idx_mol, mol in enumerate(chain(left, right)):
            lin_sys[idx_elem, idx_mol] = mol[elem]
    lin_sys[:, len(left) :] *= -1

    logger.info("Linear system of equations matrix:\n%s", lin_sys)
    nullspace: list[Matrix] = Matrix(lin_sys).nullspace(simplify=True)
    assert all(null_basis.shape[1] == 1 for null_basis in nullspace), (
        "Kernel basis should consist of column vectors."
    )
    kernel: list[list[Rational]] = [null_basis.flat() for null_basis in nullspace]

    logger.info("Nullity = %d", len(kernel))
    solutions: list[tuple[IntArr, IntArr]] = []
    for ker in kernel:
        coefs = _scale_to_integers(ker)
        logger.info("Kernel basis vector %s scaled to %s", ker, coefs)
        if np.any(coefs < 0) and np.any(coefs > 0):
            continue
        coefs = np.abs(coefs)
        cut = len(left)
        solutions.append((coefs[:cut], coefs[cut:]))
    return solutions
