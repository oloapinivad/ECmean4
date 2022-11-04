CDO ECmean4
===========

An older version of the ECmean scripts based on CDO calls is still available, and it can be found in the ``CDO`` folder. 
Here you can find the instruction to install it and run it. 

Install in a python virtual environment
---------------------------------------

It is recommended to work into a Python virtual environment `venv`, which can be created with:

.. code-block:: shell

    > python -m venv .ECmean4

The ``.ECmean4`` is just an arbitrary name, with the leading dot meaning that it is a hidden folder so that it does not mess up your $HOME.
Then you can activate the environment with:

.. code-block:: shell

    > source .ECmean4/bin/activate

Could be convenient to create an alias to activate the enviroment rapidly when required, adding in the ``.bash_profile.sh`` :: 

    alias pyECmean4='source .ECmean4/bin/activate'


Requirements
------------

The required packages are listed in ``docs/requirements.txt`` and reported for completeness here below. 
They can be easily installed with pip. ::

    cdo==1.5.6
    MetPy==1.3.0
    numpy==1.22.3
    PyYAML==6.0
    tabulate==0.8.9

`NetCDF4` is also required to accessing the files although is not formally defined as a requirement.

ECmean4 is based on two independent Python scripts which builts on the same set of functions and classes.
They both produces takes care of all the computation and produces a .txt table:

.. warning::
  The CDO version is no longer developed, and the documentation is not updated. 

global_mean.py
--------------

Main concepts
^^^^^^^^^^^^^

The ``global_mean.py`` script computes the global averages for many dynamical and physical fields. It compares the output against a set of climatological values defined in ``gm_reference.yml``, including temperature, salinity, precipitation (over land and sea), radiative fluxes and other quantities useful for model tuning.

Usage
^^^^^

Running the global mean evaluation is rather simple ::

        ./global_mean.py EXP Y1 Y2

Extra options are listed here below:

positional arguments:
  EXP                   experiment ID
  Y1                    starting year
  Y2                    final year

optional arguments:
  -h, --help            show this help message and exit
  -s, --silent          do not print anything to std output
  -t, --trend           compute trends
  -l, --line            appends also single line to a table
  -o FILE, --output FILE
                        path of output one-line table
  -m MODEL, --model MODEL
                        model name
  -v LOGLEVEL, --loglevel LOGLEVEL
                        define the level of logging.
  -j NUMPROC            number of processors to use
  -e ENSEMBLE, --ensemble ENSEMBLE
                        variant label (ripf number for cmor)
  -d, --debug           activate cdo debugging

performance_indices.py
----------------------

Main concepts
^^^^^^^^^^^^^

The ``performance_indices.py`` script computes the `Reichler and Kim Performance Indices <https://journals.ametsoc.org/view/journals/bams/89/3/bams-89-3-303.xml>`_. A few minor modification has been implemented in the original ECmean configuration so that it still work on a regular 2x2 grid (instead of using observation grid). and it is based on old climatological assessment present in the original ECmean. The climatology is defined by ``pi_climatology.yml`` script. Climatology for PI is still VERY outdated and it is being updated. 


Usage
^^^^^

Running the performance indices evaluation is rather simple ::

        ./performance_indices.py EXP Y1 Y2

Extra options are listed here below:

positional arguments:
  EXP                   experiment ID
  Y1                    starting year
  Y2                    final year

optional arguments:
  -h, --help            show this help message and exit
  -s, --silent          do not print anything to std output
  -v LOGLEVEL, --loglevel LOGLEVEL
                        define the level of logging. default: error
  -j NUMPROC            number of processors to use
  -m MODEL, --model MODEL
                        model name
  -e ENSEMBLE, --ensemble ENSEMBLE
                        variant label (ripf number for cmor)
  -k CLIMATOLOGY        which climatology you want to use (EC22 or RK08)
  -r RESOLUTION         only EC22: resolution of the climatology (r180x90 or r360x180)
  -d, --debug           activate cdo debugging

Climatology
^^^^^^^^^^^

The performance indices built on the comparison between model data and a previously computed climatology of several variables.
The original ECmean climatology is the defined as RK08, a new one being developed with newer high-resolution data is defined as EC22. 
The default climatology is the RK08, but this will be replaced in the future. 

.. warning::
	A bug is known for sea surface salinity variance, as described on the correspondent `Github Issue <https://github.com/oloapinivad/ECmean4/issues/8>`_ Please be aware the this PI is affected. 



