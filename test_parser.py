"""
Unit tests for chemical parser.

Copyright 2022. Andrew Wang.
"""
from typing import Dict
from parse import parse


def test_element():
    """Unit tests for elemental compounds."""
    check('I', {'I': 1})
    check('Os', {'Os': 1})
    check('N3', {'N': 3})
    check('Cl2', {'Cl': 2})


def test_compound():
    """Unit tests for simple compounds."""
    check('KCl', {'K': 1, 'Cl': 1})
    check('NO', {'N': 1, 'O': 1})
    check('Na2SO4', {'Na': 2, 'S': 1, 'O': 4})
    check('Na2SO4', {'Na': 2, 'S': 1, 'O': 4})
    check('NaHCO3', {'Na': 1, 'H': 1, 'C': 1, 'O': 3})
    check('C6H12O6', {'C': 6, 'H': 12, 'O': 6})
    check('MoSi2', {'Mo': 1, 'Si': 2})


def test_paren():
    """Unit tests with non-nested parenthesis."""
    check('PtCl2(NH3)2', {'Pt': 1, 'Cl': 2, 'N': 2, 'H': 6})
    check('Fe2(SO4)3', {'Fe': 2, 'S': 3, 'O': 12})
    check('Fe(H2O)4(OH)2', {'Fe': 1, 'H': 10, 'O': 6})
    check('(Li3)(NaO2)', {'Li': 3, 'Na': 1, 'O': 2})
    check('(CH3)(CH2)2', {'C': 3, 'H': 7})
    check('Ca3(PO4)2Hg(NH3)', {'Ca': 3, 'P': 2,
          'O': 8, 'Hg': 1, 'N': 1, 'H': 3})


def test_nested():
    """Unit tests with nested parenthesis."""
    check('TiCl2((CH3)2PCH2CH2P(CH3)2)2', {
          'Ti': 1, 'Cl': 2, 'C': 12, 'P': 4, 'H': 32})
    check('NaRb5(PuO4(OH)2)2', {'Na': 1, 'Rb': 5, 'Pu': 2, 'O': 12, 'H': 4})
    check('(PuO4(OH)2)NaM2', {'Pu': 1, 'O': 6, 'H': 2, 'Na': 1, 'M': 2})
    check('CH2(CH2CH(C6H5))2CH2', {'C': 18, 'H': 20})


def check(mol: str, expected: Dict[str, int]):
    """Assert that the result of the parser is equivalent to expected."""
    actual = parse(mol)
    assert dict(actual) == expected
