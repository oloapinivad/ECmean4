Xarray Scripts
==============
ECmean4 is based on two independent Python scripts which builts on the same set of functions and classes ``xr_ecmean.py``.
They both produces takes care of all the computation and produces a .txt table:

xr_global_mean.py
-----------------

Main concepts
^^^^^^^^^^^^^

The ``xr_global_mean.py`` script computes the global averages for many dynamical and physical fields. It compares the output against a set of climatological values defined in ``gm_reference.yml``, including temperature, salinity, precipitation (over land and sea), radiative fluxes and other quantities useful for model tuning.

Usage
^^^^^

Running the global mean evaluation is rather simple ::

        ./xr_global_mean.py EXP Y1 Y2

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

xr_performance_indices.py
-------------------------

Main concepts
^^^^^^^^^^^^^

The ``xr_performance_indices.py`` script computes the `Reichler and Kim Performance Indices <https://journals.ametsoc.org/view/journals/bams/89/3/bams-89-3-303.xml>`_. A few minor modification has been implemented in the original ECmean configuration so that it still work on a regular 2x2 grid (instead of using observation grid). and it is based on old climatological assessment present in the original ECmean. The climatology is defined by ``pi_climatology.yml`` script. Climatology for PI is still VERY outdated and it is being updated. 


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
  -d, --debug           activate cdo debugging
  -k CLIMATOLOGY        which climatology you want to use (EC22 or RK08)
  -r RESOLUTION         only EC22: resolution of the climatology (r180x90 or r360x180)


Climatology
^^^^^^^^^^^

The performance indices built on the comparison between model data and a previously computed climatology of several variables.
The original ECmean climatology is the defined as RK08, a new one being developed with newer high-resolution data is defined as EC22. 
The default climatology is the RK08, but this will be replaced in the future. 

.. warning::
	A bug is known for sea surface salinity variance, as described on the correspondent `Github Issue <https://github.com/oloapinivad/ECmean4/issues/8>`_ Please be aware the this PI is affected. 

CMOR compatibility
------------------

It is possible to use these tools also to analyze CMOR files for CMIP5 or CMIP6. This assumes a standard ESGF directory structure but you can change it by modifying the corresponding interface files ``interfaces/interface_CMIP6.yml`` and ``interfaces/interface_CMIP6.yml``.
In order to allow masking you will need the `sftlf`, `sftof` and `areacello` variables for you experiment of interest too.

Usage example for CMIP5::

        ./xr_global_mean.py historical 1990 1999 -j 12 -m EC-EARTH -e r1i1p1 -i CMIP5

will compute performance indices for member r1i1p1 of the EC-EARTH model in the CMIP5 historical experiment.

Usage example for CMIP6::

        ./xr_performance_indices.py historical 1990 1999 -j 12 -m EC-Earth3 -e r1i1p1f1 -i CMIP6 -k EC22 -r r360x180

The same as above, but for the CMIP6 EC-Earth3 model. In this case the comparison is with the newer EC22 climatology at high r360x180 resolution.