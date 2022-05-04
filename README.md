[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
![Basic Pytest](https://github.com/github/docs/actions/workflows/pytest.yml/badge.svg)

![ECmean4](docs/ecmean4_smallest.png)

*An EC-Earth4 basic evaluation tool*

This is the painful road to switch from bash to python3 as a scripting interface.
As in the previous ECmean, all computations and weights are done with CDO. However, python3 is taking care of the script in place of bash, and YML files are using to handle all the configuration/metadata required. 

## Scripts

Two scripts takes care of all the computation and produces a txt table:

- **global_mean.py**: computes the global averages for many dynamical and physical fields. It compares the output against a set of climatological values defined in `reference.yml`
- **performance_indices.py**: computes Reichler and Kim Performance Indices. It work on a regular 2x2 grid and it is based on old climatological assessment present in the original ECmean. The climatology is defined by `pi_climatology.yml` script. Climatology for PI is still VERY outdated.

## Configuration

A template configuration file is included in the repo, including folder definition and all the details that can be manually adjusted
You will need to copy `cp config.tmpl config.yml` and then edit the required folder before runninn
g the script`

Conversion from model variables - as well as file structure - is handled by a model depedent file `interface_ece4.yml` that can be potentially updated to other climate models.

## Python dependencies and installation

Evident dependencies are python3 packages cdo, PyYAML, netCDF, metpy (pint) and tabulate. 
Python 3.7+ is required. 

Recommend to work into a python virtual environment

`python3 -m venv .ECmean4`

Then install everything with pip

## Extra

More will come, sooner or later.

