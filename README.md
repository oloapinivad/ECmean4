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

More will come, sooner or later

