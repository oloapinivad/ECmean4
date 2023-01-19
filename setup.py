#!/usr/bin/env python

from distutils.core import setup

setup(name='ECmean4',
      version='0.0.1',
      description='ECmean4 lightweight evaluation tool',
      author='Paolo Davini, Jost bon Hardenberg',
      author_email='p.davini@isac.cnr.it',
      url='https://github.com/oloapinivad/ECmean4',
      packages=['ecmean'],
      entry_points={
            "console_scripts": [
                "global_mean = ecmean.global_mean:gm_entry_point",
            ],
            }
     )

