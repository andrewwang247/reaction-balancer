"""
Chemical reaction balance solver.

Copyright 2022. Andrew Wang.
"""
from typing import DefaultDict, Iterable, List, Tuple
from itertools import chain
import numpy as np
from sympy import Matrix


def __distinct_elems(mols: List[DefaultDict[str, int]]) -> List[str]:
    """Get the distinct elements that form the molecules."""
    elems = set()
    for mol in mols:
        for key in mol:
            elems.add(key)
    return list(elems)


def balance(left: List[DefaultDict[str, int]],
            right: List[DefaultDict[str, int]],
            verbose: bool) \
        -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """Balance the parsed left and ride sides."""
    elems = __distinct_elems(left + right)
    if verbose:
        print('Distinct elements:', ', '.join(elems))
    lin_sys = np.zeros((len(elems), len(left) + len(right)), dtype=int)
    for idx_elem, elem in enumerate(elems):
        for idx_mol, mol in enumerate(chain(left, right)):
            lin_sys[idx_elem, idx_mol] = mol[elem]
    lin_sys[:, len(left):] *= -1
    if verbose:
        print('Linear system of equations matrix', lin_sys, sep='\n')
    kernel = Matrix(lin_sys).nullspace(simplify=True)
    if verbose:
        print('Nullity =', len(kernel))
    for ker in kernel:
        numers = np.array([num.as_numer_denom()[0] for num in ker])
        denoms = np.array([num.as_numer_denom()[1] for num in ker])
        lcm = np.lcm.reduce(denoms)
        coefs = numers * lcm / denoms
        coefs /= np.gcd.reduce(coefs)
        if np.any(coefs < 0) and np.any(coefs > 0):
            continue
        coefs = np.abs(coefs)
        if verbose:
            print('Kernel basis', ker, 'simplified to', coefs)
        yield coefs[:len(left)], coefs[len(left):]
