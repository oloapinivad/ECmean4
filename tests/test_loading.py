"""Tests for data loading and work split (#155)"""

import ecmean.libs.files
from ecmean import GlobalMean
from ecmean.libs.files import load_variable
from ecmean.libs.general import runtime_weights, weight_split


def test_weight_split_balanced():
    """Variables are assigned to the least loaded worker (8/10 with the previous lexicographic split)"""
    variables = ["ta", "pme", "net_sfc", "net_toa", "tas"]
    chunks = weight_split(variables, 2)
    weights = runtime_weights(variables)
    assert sorted(sum(weights[v] for v in chunk) for chunk in chunks) == [9, 9]
    assert sorted(v for chunk in chunks for v in chunk) == sorted(variables)


def test_load_variable_opens_files_once(monkeypatch):
    """Input files are opened only once per variable"""
    calls = []
    original = ecmean.libs.files.open_input_files

    def counting(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)

    monkeypatch.setattr(ecmean.libs.files, "open_input_files", counting)
    gm = GlobalMean("cpld", 1990, 1990, config="tests/config.yml", reference="EC26-HIST")
    gm.prepare()
    field, units = load_variable("tas", gm.face, gm.diag)
    assert len(calls) == 1
    assert field.notnull().any()
    assert units == "K"
