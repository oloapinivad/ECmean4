Performance Indices
===================

Main concepts
^^^^^^^^^^^^^

The ``performance_indices.py`` script computes the `Reichler and Kim Performance Indices <https://journals.ametsoc.org/view/journals/bams/89/3/bams-89-3-303.xml>`_. 
Some minor differences from the original definition has been introduced, so that the PIs are computed on a common grid rather than on the original grid.
From the original definition a few improvements has been introduced, producing the PIs also for a set of selected regions and seasons. 

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
  -c config.yml, --config config.yml
                        set up a specific configuration file (config.yml is default)
  -i INTERFACE, --interface INTERFACE
                        set up a specific interface file (override config.yml)
  -m MODEL, --model MODEL
                        model name
  -e ENSEMBLE, --ensemble ENSEMBLE
                        variant label (ripf number for cmor)
  -d, --debug           activate cdo debugging
  -k CLIMATOLOGY        which climatology you want to use (EC23: default, RK08 alternative under update)
  -r RESOLUTION         only EC22: resolution of the climatology (r180x90 or r360x180)

Example 
^^^^^^^

Usage example for CMIP6::

        ./performance_indices.py historical 1990 1999 -j 12 -m EC-Earth3 -e r1i1p1f1 -i CMIP6 -k EC22 -r r360x180

The same as above, but for the CMIP6 EC-Earth3 model. In this case the comparison is with the newer EC22 climatology at high r360x180 resolution.


Output
^^^^^^

The result is produced in a form a YAML file that can be stored for later comparison. 
Most importantly, a figure is produced showing a score card for the different regions, variables and seasons.
This is computed as the ration between the model PI and the average value estimated over the ensemble of CMIP6 models. 
An example of the the output for a single year of the EC-Earth3 historical simulation is shown here below.

.. figure:: pitestfigure.png
   :align: center
   :width: 600px
   :alt: PI for ECearth3

   An example for a single year of the EC-Earth3 historical r1i1p1f1 simulation.

Such plot is currently available for the EC23 climatology only, which is currently computed on a 30-year time window from 1990 to 2019.
Similarly, season-dependent computation are available only for EC23.
Details on the field used are reported in the `climatology/EC23/pi_climatology_EC23.yml` file.


Climatology
^^^^^^^^^^^

The performance indices built on the comparison between model data and a previously computed climatology of several variables.
The original ECmean climatology was the defined as RK08. However, a new has been developed with high-resolution data and is now defined as EC23. 
The default climatology is the EC23. An intermadiate versio knowns as EC22 is available but not recommended and will be removed soon.

.. warning::
	A bug is known for sea surface salinity variance, as described on the correspondent `Github Issue <https://github.com/oloapinivad/ECmean4/issues/8>`_ Please be aware the this PI is affected. 

Climatology is computed by the `py-climatology-create.py` script, which is included in the repository for documentation.
It is based on a basic YAML file which is based on the machine where the climatology has been developed and it is used to identifiy all the required data folder and names. 
The tool loop over the variable and produces the yearly and seasonal average of the climate, as well as the interannual variance. 
To avoid that grid points with irrealistic low variance affect the computation of the PIs, a filter based on the log10 5 sigma is introduced.

Once the climatology is created, the script `cmip6-clim-evaluate.py` is used to run iteratively on a set of 10 CMIP6 models and later to compute the multi model mean of the PIs (for each region and season).
This is later used to provide a ratio between the original PI and the CMIP6 ensemble. 

