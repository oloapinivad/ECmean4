Global Mean
===========

Main concepts
^^^^^^^^^^^^^

The ``global_mean.py`` script computes the global averages for many dynamical and physical fields. 
It compares the output against a set of climatological values defined in ``gm_reference.yml``, including temperature, salinity, precipitation (integrals over land and sea), radiative fluxes and other quantities useful for model tuning.

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
  -t, --trend           compute trends on multiple years
  -l, --line            appends also single line to a table
  -c config.yml, --config       
                        config.yml set up a specific configuration file (config.yml is default)
  -i INTERFACE, --interface INTERFACE   set up a specific interface file (override config.yml)
  -o FILE, --output FILE        path of output one-line table
  -m MODEL, --model MODEL       model name (override config.yml)
  -v LOGLEVEL, --loglevel LOGLEVEL      define the level of logging. Default is warning.
  -j NUMPROC                    number of processors to use
  -e ENSEMBLE, --ensemble ENSEMBLE      variant label (ripf number for cmor)


Output
^^^^^^

A txt table including the comparison with some predefined dataset.
