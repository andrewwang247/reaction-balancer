"""Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""

from typing import TYPE_CHECKING

from src import solve

if TYPE_CHECKING:
    from .conftest import Equation


def test_balance(equation: Equation) -> None:
    """Assert that the balanced equations are correct."""
    solutions = list(solve(equation.left_mols, equation.right_mols))
    assert len(solutions) == 1, "Solution should be unique."
    left_coefs, right_coefs = solutions[0]
    assert left_coefs.tolist() == equation.left_coefs
    assert right_coefs.tolist() == equation.right_coefs
