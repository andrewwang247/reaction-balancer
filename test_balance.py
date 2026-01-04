"""
Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""
from typing import Dict, List, Union
from json import load
from pytest import mark
from balance import solve


def _get_test_eqns(
        filename: str) -> List[Dict[str, Union[List[str], List[int]]]]:
    """Parse the test cases JSON."""
    with open(filename, encoding='UTF-8') as fp:
        return load(fp)


@mark.parametrize('equation', _get_test_eqns('tst/equations.json'))
def test_balance(equation: Dict[str, Union[List[str], List[int]]]):
    """Assert that the balanced equations are correct."""
    left_mols = equation['left_mols']
    right_mols = equation['right_mols']
    assert all(isinstance(mol, str) for mol in left_mols)
    assert all(isinstance(mol, str) for mol in right_mols)
    left_coefs_expected = equation['left_coefs']
    right_coefs_expected = equation['right_coefs']
    assert all(isinstance(coef, int) for coef in left_coefs_expected)
    assert all(isinstance(coef, int) for coef in right_coefs_expected)
    solutions = list(solve(left_mols, right_mols, False))  # type:  ignore
    assert len(solutions) == 1, 'Solution should be unique.'
    left_coefs_actual, right_coefs_actual = solutions[0]
    assert all(left_coefs_actual == left_coefs_expected)
    assert all(right_coefs_actual == right_coefs_expected)
