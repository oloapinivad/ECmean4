"""Tests for the formula tokenizer and evaluator"""

import numpy as np
import pytest
import xarray as xr

from ecmean.libs.files import _get_variables_to_load
from ecmean.libs.formula import _eval_formula, formula_variables

DATASET = xr.Dataset({"a": ("x", np.array([4.0, 8.0])), "b": ("x", np.array([1.0, 2.0])), "c": ("x", np.array([2.0, 2.0]))})


@pytest.mark.parametrize(
    "formula, expected",
    [
        ("a", [4.0, 8.0]),
        ("a+b", [5.0, 10.0]),
        ("a - b", [3.0, 6.0]),
        ("a*0.5", [2.0, 4.0]),
        ("a*1e3", [4000.0, 8000.0]),
        ("a*2.5E-1", [1.0, 2.0]),
        ("(a-b)*2", [6.0, 12.0]),
        ("a-b-c", [1.0, 4.0]),
        ("a/b/c", [2.0, 2.0]),
        ("a-(b-c)", [5.0, 8.0]),
        ("a+b*c", [6.0, 12.0]),
        ("-a+b", [-3.0, -6.0]),
        ("2*(a+b)/c", [5.0, 10.0]),
    ],
)
def test_eval_formula(formula, expected):
    """Check precedence, parentheses, spaces and numeric literals"""
    np.testing.assert_allclose(_eval_formula(formula, DATASET).values, expected)


def test_formula_variables():
    """Variables are extracted once, ignoring numbers and exponents"""
    assert formula_variables("(a - b)*2 + a*1e3") == ["a", "b"]
    assert _get_variables_to_load("x", {"variables": {"x": {"derived": "a*0.5+b"}}}) == ["a", "b"]


@pytest.mark.parametrize("formula", ["a$b", "a+", "(a+b", "a+b)", "a b", "", "*a"])
def test_unsupported_formula(formula):
    """Unsupported or malformed syntax gives a clear ValueError"""
    with pytest.raises(ValueError):
        _eval_formula(formula, DATASET)
