"""Configure pytest fixtures.

Copyright 2026. Andrew Wang.
"""

from dataclasses import dataclass
from json import load
from pathlib import Path
from typing import Any

import pytest

from src import Elements, parse

_RESOURCE_DIR = Path("resources")


@dataclass
class Molecule:
    """JSON structure for molecules."""

    name: str
    elements: dict[str, int]


@dataclass
class Equation:
    """JSON structure for equations."""

    left_mols: list[Elements]
    right_mols: list[Elements]
    left_coefs: list[int]
    right_coefs: list[int]


def _get_resource(resource: str) -> list[dict[str, Any]]:
    """Retrieve a resource based on its name."""
    fpath = _RESOURCE_DIR / f"{resource}.json"
    assert fpath.is_file()
    with fpath.open(encoding="UTF-8") as fp:
        rs: list[dict[str, Any]] = load(fp)
    return rs


@pytest.fixture(scope="session", params=_get_resource("molecules"))
def molecule(request: pytest.FixtureRequest) -> Molecule:
    """Configure molecule data."""
    rs: dict[str, Any] = request.param
    mol_name: str = rs["molecule"]
    elements: dict[str, int] = rs["elements"]
    return Molecule(mol_name, elements)


@pytest.fixture(scope="session", params=_get_resource("equations"))
def equation(request: pytest.FixtureRequest) -> Equation:
    """Configure equation data."""
    rs: dict[str, Any] = request.param

    left_mols = [parse(mol) for mol in rs["left_mols"]]
    right_mols = [parse(mol) for mol in rs["right_mols"]]
    left_coefs: list[int] = rs["left_coefs"]
    right_coefs: list[int] = rs["right_coefs"]
    return Equation(left_mols, right_mols, left_coefs, right_coefs)
