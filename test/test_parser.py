"""
Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""
from typing import cast, Dict, List, TypedDict
from json import load
from pytest import mark
from src import parse


class Molecule(TypedDict):
    """JSON structure for molecules."""
    molecule: str
    elements: Dict[str, int]


def _get_test_molecules() -> List[Molecule]:
    with open('test/molecules.json', encoding='UTF-8') as fp:
        return cast(List[Molecule], load(fp))


@mark.parametrize('case', _get_test_molecules())
def test_parser(case: Molecule) -> None:
    """Assert that the result of the parser is equivalent to expected."""
    molecule = case['molecule']
    expected_elements = case['elements']
    actual_elements = parse(molecule)
    assert dict(actual_elements) == expected_elements
