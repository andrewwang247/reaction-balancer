"""
Chemical reaction balance solver.

Copyright 2026. Andrew Wang.
"""
from typing import DefaultDict, Iterable, List, Tuple
from itertools import chain
import numpy as np
from sympy import Matrix, Rational  # type: ignore
from parse import parse


def __distinct_elems(mols: List[DefaultDict[str, int]]) -> List[str]:
    """Get the distinct elements that form the molecules."""
    elems = set()
    for mol in mols:
        for key in mol:
            elems.add(key)
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
            right: List[DefaultDict[str, int]],
            verbose: bool) \
        -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """Balance the parsed left and ride sides."""
    elems = __distinct_elems(left + right)
    if verbose:
        print(f'Distinct elements ({len(elems)}):', elems)
    lin_sys = np.zeros((len(elems), len(left) + len(right)), dtype=int)
    for idx_elem, elem in enumerate(elems):
        for idx_mol, mol in enumerate(chain(left, right)):
            lin_sys[idx_elem, idx_mol] = mol[elem]
    lin_sys[:, len(left):] *= -1
    if verbose:
        print('Linear system of equations matrix:', lin_sys, sep='\n')
    nullspace: List[Matrix] = Matrix(lin_sys).nullspace(simplify=True)
    assert all(null_basis.shape[1] == 1 for null_basis in nullspace), \
        'Kernel basis should consist of column vectors.'
    kernel: List[List[Rational]] = [null_basis.flat()
                                    for null_basis in nullspace]
    if verbose:
        print(f'Nullity = {len(kernel)}. Kernel basis:')
        for ker in kernel:
            print(f'\t{ker}')
    for ker in kernel:
        coefs = __scale_to_integers(ker)
        if np.any(coefs < 0) and np.any(coefs > 0):
            continue
        coefs = np.abs(coefs)
        yield coefs[:len(left)], coefs[len(left):]


def solve(left: Iterable[str], right: Iterable[str],
          verbose: bool) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """Balance left and right sides of chemical equation."""
    left_mols = [parse(mol) for mol in left]
    right_mols = [parse(mol) for mol in right]
    if verbose:
        print('Molecules (L):')
        for mol, elem_counts in zip(left, left_mols):
            print(f'\t{mol}:', dict(elem_counts))
        print('Molecules (R):')
        for mol, elem_counts in zip(right, right_mols):
            print(f'\t{mol}:', dict(elem_counts))
    return balance(left_mols, right_mols, verbose)
