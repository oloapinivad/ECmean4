[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/oloapinivad/ECmean4/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Documentation Status](https://readthedocs.org/projects/ecmean4/badge/?version=latest)](https://ecmean4.readthedocs.io/en/latest/?badge=latest)
[![Basic PyTest](https://github.com/oloapinivad/ECmean4/actions/workflows/pytest.yml/badge.svg)](https://github.com/oloapinivad/ECmean4/actions/workflows/pytest.yml)


![ECmean4](docs/ecmean4_smallest.png)

*An EC-Earth4 basic evaluation tool*

ECmean4 is a lightweight parallelized tool for evaluation of basic properties of Global Climate Models, such as global mean, pattern correlation and climate model performance indices.

It builts on the original ECmean which has been used for EC-Earth2 and EC-Earth3 evaluation, but it uses Python3 as a scripting language to perform lazy calls to CDO and makes use of YML configuration files. It works both on raw EC-Earth4 output and on CMOR model output from CMIP5 and CMIP6.

A complete [ReadTheDocs Documentation](https://ecmean4.readthedocs.io/en/latest/index.html) is being developed, please refer to this.
