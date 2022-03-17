[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

# ECmean4
EC-Earth basic evaluation tool

This is the painful road to switch from bash to python as a scripting interface.
As in the previous ECmean, all the inner computations and weights are done with CDO

Two scripts takes care of all the computation, one for global mean and the other for the Reichler and Kim Performance Indices

Interface to configuration file or climatoligical assessment is provided by YAML files
Climatology for PI is still VERY outdated.

Evident dependencies are python3 packages cdo, PyYAML, netCDF4 and tabulate

Recommend to work into a python virtual environment

`python3 -m venv .ECmean4`

Then install everything with pip

A template configuration file is included in the repo, including folder definition and all the details that can be manually adjusted
You will need to copy `cp config.tmpl config.yml` and then edit the required folder before running the script`

More will come, sooner or later

