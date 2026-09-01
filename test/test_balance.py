"""Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""

from json import load
from pathlib import Path
from typing import TypedDict, cast

import pytest

from src import solve


class Equation(TypedDict):
    """JSON structure for equations."""

    left_mols: list[str]
    right_mols: list[str]
    left_coefs: list[int]
    right_coefs: list[int]


def _get_test_equations() -> list[Equation]:
    """Parse the test cases JSON."""
    with Path("test/equations.json").open(encoding="UTF-8") as fp:
        return cast("list[Equation]", load(fp))


@pytest.mark.parametrize("equation", _get_test_equations())
def test_balance(equation: Equation) -> None:
    """Assert that the balanced equations are correct."""
    left_mols = equation["left_mols"]
    right_mols = equation["right_mols"]
    left_coefs_expected = equation["left_coefs"]
    right_coefs_expected = equation["right_coefs"]
    solutions = list(solve(left_mols, right_mols))
    assert len(solutions) == 1, "Solution should be unique."
    left_coefs_actual, right_coefs_actual = solutions[0]
    assert left_coefs_actual.tolist() == left_coefs_expected
    assert right_coefs_actual.tolist() == right_coefs_expected
