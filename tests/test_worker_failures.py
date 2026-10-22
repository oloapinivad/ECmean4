"""Tests for robustness of the workers: a failing variable must not drop the others (#151)"""

import importlib
import math
import os

import pytest
import yaml

from ecmean import GlobalMean, PerformanceIndices

# the package exports functions with the same name of the modules, so import the modules explicitly
gm_module = importlib.import_module("ecmean.global_mean")
pi_module = importlib.import_module("ecmean.performance_indices")

VARS = ["tas", "psl"]


def _failing_formula(module, monkeypatch, failvar="psl"):
    """Make formula_wrapper raise for a single variable, workers are forked so they inherit the patch"""
    original = module.formula_wrapper

    def wrapper(var, face, xfield):
        if var == failvar:
            raise RuntimeError(f"Simulated failure for {var}")
        return original(var, face, xfield)

    monkeypatch.setattr(module, "formula_wrapper", wrapper)


def test_global_mean_failing_variable(monkeypatch, tmp_path):
    """A failing variable is stored as NaN, the other variables of the same worker survive"""
    _failing_formula(gm_module, monkeypatch)
    gm = GlobalMean("cpld", 1990, 1990, config="tests/config.yml", numproc=1, reference="EC26-HIST")
    gm.diag.var_all = VARS
    gm.prepare()
    gm.run()
    yamlfile = tmp_path / "gm.yml"
    gm.store(yamlfile=yamlfile, tablefile=tmp_path / "gm.txt")

    with open(yamlfile, "r", encoding="utf8") as file:
        data = yaml.safe_load(file)
    assert not math.isnan(data["tas"]["ALL"]["Global"])
    assert math.isnan(data["psl"]["ALL"]["Global"])


def test_performance_indices_failing_variable(monkeypatch, tmp_path):
    """A failing variable is stored as NaN, the other variables of the same worker survive"""
    _failing_formula(pi_module, monkeypatch)
    pi = PerformanceIndices("cpld", 1990, 1990, config="tests/config.yml", numproc=1, climatology="EC26-HIST")
    pi.diag.field_all = VARS
    pi.prepare()
    pi.run()
    yamlfile = tmp_path / "pi.yml"
    pi.store(yamlfile=yamlfile)

    with open(yamlfile, "r", encoding="utf8") as file:
        data = yaml.safe_load(file)
    assert not math.isnan(data["tas"]["ALL"]["Global"])
    assert math.isnan(data["psl"]["ALL"]["Global"])


def test_global_mean_dead_worker(monkeypatch):
    """A worker that dies without raising (e.g. killed by the OS) makes run() fail explicitly"""

    def die(var, face, xfield):
        os._exit(1)

    monkeypatch.setattr(gm_module, "formula_wrapper", die)
    gm = GlobalMean("cpld", 1990, 1990, config="tests/config.yml", numproc=1, reference="EC26-HIST")
    gm.diag.var_all = VARS
    gm.prepare()
    with pytest.raises(RuntimeError, match="tas"):
        gm.run()
