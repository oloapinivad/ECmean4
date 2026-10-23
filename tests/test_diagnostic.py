"""Unit tests for the configuration checks of ecmean.libs.diagnostic."""

import copy

import pytest
import yaml

from ecmean.libs.diagnostic import Diagnostic

BASE_CONFIG = {
    "interface": "EC-Earth4",
    "model": {"name": "EC-Earth4"},
    "dirs": {"exp": "tests/data", "tab": "tests/table", "fig": "tests/figures"},
    "global_mean": {
        "variables": {"atm": ["tas"]},
        "regions": ["Global", "Tropical"],
        "seasons": ["ALL", "JJA"],
        "reference": "EC26-CMIP",
    },
    "performance_indices": {
        "variables": {"atm2d": ["tas"]},
        "regions": ["Global", "Tropical"],
        "seasons": ["ALL", "JJA"],
        "climatology": "EC26-CMIP",
    },
}


def _diagnostic(funcname, config=None, **kwargs):
    return Diagnostic("exp", 1990, 1990, BASE_CONFIG if config is None else config, funcname, **kwargs)


@pytest.mark.parametrize(
    "key, value, message",
    [("regions", ["Tropical"], "'Global' region"), ("seasons", ["JJA"], "'ALL' season")],
)
def test_global_all_mandatory_global_mean(key, value, message):
    config = copy.deepcopy(BASE_CONFIG)
    config["global_mean"][key] = value
    with pytest.raises(ValueError, match=message):
        _diagnostic("GlobalMean", config)


def test_global_all_not_required_for_performance_indices():
    config = copy.deepcopy(BASE_CONFIG)
    config["performance_indices"]["regions"] = ["NH", "SH"]
    config["performance_indices"]["seasons"] = ["JJA"]
    assert _diagnostic("PerformanceIndices", config).regions == ["NH", "SH"]


def test_versions_are_separate():
    diag = _diagnostic("GlobalMean")
    assert diag.version == "*"
    assert diag.ecmean_version != "*"


def test_resolution_default():
    assert _diagnostic("PerformanceIndices").resolution == "r360x180"


def test_resolution_from_config():
    config = copy.deepcopy(BASE_CONFIG)
    config["performance_indices"]["resolution"] = "r360x180"
    assert _diagnostic("PerformanceIndices", config).resolution == "r360x180"


@pytest.mark.parametrize("from_argument", [False, True])
def test_resolution_not_supported(from_argument):
    config = copy.deepcopy(BASE_CONFIG)
    kwargs = {}
    if from_argument:
        kwargs["resolution"] = "r90x45"
    else:
        config["performance_indices"]["resolution"] = "r90x45"
    with pytest.raises(ValueError, match="not supported"):
        _diagnostic("PerformanceIndices", config, **kwargs)


def test_default_config_missing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match="config.yml"):
        _diagnostic("GlobalMean", config="")


def test_default_config_in_current_directory(tmp_path, monkeypatch):
    (tmp_path / "config.yml").write_text(yaml.safe_dump(BASE_CONFIG))
    monkeypatch.chdir(tmp_path)
    assert _diagnostic("GlobalMean", config="").modelname == "EC-Earth4"
