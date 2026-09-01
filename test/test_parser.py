"""Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""

from json import load
from pathlib import Path
from typing import TypedDict, cast

import pytest

from src import parse


class Molecule(TypedDict):
    """JSON structure for molecules."""

    molecule: str
    elements: dict[str, int]


def _get_test_molecules() -> list[Molecule]:
    with Path("test/molecules.json").open(encoding="UTF-8") as fp:
        return cast("list[Molecule]", load(fp))


@pytest.mark.parametrize("case", _get_test_molecules())
def test_parser(case: Molecule) -> None:
    """Assert that the result of the parser is equivalent to expected."""
    molecule = case["molecule"]
    expected_elements = case["elements"]
    actual_elements = parse(molecule)
    assert dict(actual_elements) == expected_elements
