"""
Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""
from typing import cast, List, TypedDict
from json import load
from pytest import mark
from src import solve


class Equation(TypedDict):
    """JSON structure for equations."""
    left_mols: List[str]
    right_mols: List[str]
    left_coefs: List[int]
    right_coefs: List[int]


def _get_test_equations() -> List[Equation]:
    """Parse the test cases JSON."""
    with open('test/equations.json', encoding='UTF-8') as fp:
        return cast(List[Equation], load(fp))


@mark.parametrize('equation', _get_test_equations())
def test_balance(equation: Equation) -> None:
    """Assert that the balanced equations are correct."""
    left_mols = equation['left_mols']
    right_mols = equation['right_mols']
    left_coefs_expected = equation['left_coefs']
    right_coefs_expected = equation['right_coefs']
    solutions = list(solve(left_mols, right_mols))
    assert len(solutions) == 1, 'Solution should be unique.'
    left_coefs_actual, right_coefs_actual = solutions[0]
    assert left_coefs_actual.tolist() == left_coefs_expected
    assert right_coefs_actual.tolist() == right_coefs_expected
