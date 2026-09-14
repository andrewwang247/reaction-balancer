"""Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""

from typing import TYPE_CHECKING

from src import parse

if TYPE_CHECKING:
    from .conftest import Molecule


def test_parser(molecule: Molecule) -> None:
    """Assert that the result of the parser is equivalent to expected."""
    actual_elements = parse(molecule.name)
    assert dict(actual_elements) == molecule.elements
