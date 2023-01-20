#!/usr/bin/env python

from distutils.core import setup

setup(name='ECmean4',
      version='0.0.1',
      description='ECmean4 Global Climate Model lightweight evaluation tool',
      author='Paolo Davini, Jost bon Hardenberg',
      author_email='p.davini@isac.cnr.it',
      url='https://github.com/oloapinivad/ECmean4',
      python_requires='>3.7, <3.11',
      packages=['ecmean'],
      install_requires=[
        'numpy',
        'xarray',
        'netcdf4',
        'dask',
        'esmpy',
        'matplotlib',
        'metpy',
        'tabulate',
        'cfgrib',
        'seaborn'
      ],
      entry_points={
            "console_scripts": [
                "global_mean = ecmean.global_mean:gm_entry_point",
                "performance_indices = ecmean.performance_indices:pi_entry_point",
            ],
            }
     )

