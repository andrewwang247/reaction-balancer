"""
Unit tests for chemical parser.

Copyright 2026. Andrew Wang.
"""
from typing import Dict, List, Union
from json import load
from pytest import mark
from parse import parse


def _get_test_mols(
        filename: str) -> List[Dict[str, Union[str, Dict[str, int]]]]:
    with open(filename, encoding='UTF-8') as fp:
        return load(fp)


@mark.parametrize('case', _get_test_mols('tst/molecules.json'))
def test_parser(case: Dict[str, Union[str, Dict[str, int]]]):
    """Assert that the result of the parser is equivalent to expected."""
    assert isinstance(case['molecule'], str)
    molecule: str = case['molecule']
    assert isinstance(case['elements'], dict)
    expected_elements: Dict[str, int] = case['elements']
    actual_elements = parse(molecule)
    assert dict(actual_elements) == expected_elements
