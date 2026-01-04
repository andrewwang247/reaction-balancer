"""
Chemical formula recursive parser.

Copyright 2026. Andrew Wang.
"""
from typing import DefaultDict
from itertools import chain
from collections import defaultdict
import re

PATT = re.compile(r'([A-Z][a-z]*)(\d*)')
SUB_PATT = re.compile(r'\d+')


def _find_closing_paren(mol: str, start_idx: int) -> int:
    """Find a closing paren given the starting idx of opening paren."""
    stack_count = 1
    for idx in range(start_idx + 1, len(mol)):
        if mol[idx] == '(':
            stack_count += 1
        elif mol[idx] == ')':
            stack_count -= 1
        if stack_count == 0:
            return idx
    raise ValueError(
        f'Could not find closing paren for {mol} at index {start_idx}.')


def parse(mol: str) -> DefaultDict[str, int]:
    """Parse the molecule into constituent element counts."""
    elements: DefaultDict[str, int] = defaultdict(int)
    left_paren = mol.find('(')
    if left_paren == -1:
        for match in PATT.finditer(mol):
            elem, sub = match.group(1, 2)
            elements[elem] += int(sub) if sub else 1
    else:
        right_paren = _find_closing_paren(mol, left_paren)
        paren = parse(mol[left_paren + 1: right_paren])
        sub_match = SUB_PATT.match(mol[right_paren + 1:])
        sub_paren = int(sub_match.group(0)) if sub_match else 1
        sub_offset = sub_match.end(0) if sub_match else 0
        left = parse(mol[0:left_paren])
        right = parse(mol[right_paren + sub_offset + 1:])
        for elem, count in chain(left.items(), right.items()):
            elements[elem] += count
        for elem, count in paren.items():
            elements[elem] += sub_paren * count
    return elements
