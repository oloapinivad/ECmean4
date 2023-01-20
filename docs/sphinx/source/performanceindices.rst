Performance Indices
===================

Main concepts
^^^^^^^^^^^^^

The ``performance_indices`` command is based on the ``performance_indices.py`` script which computes the `Reichler and Kim Performance Indices <https://journals.ametsoc.org/view/journals/bams/89/3/bams-89-3-303.xml>`_, usually known as PIs. 
Some minor differences from the original definition has been introduced, so that the PIs are computed on a common (user defined) grid rather than on the original grid.
From the original definition a few improvements has been introduced, producing the PIs also for a set of selected regions and seasons. 

PIs are the root mean square error of a 2D field normalized by the interannual variance estimated from the observations. Larger values implies worse performance of the climate models.

Usage
^^^^^

Running the performance indices evaluation is rather simple ::

        performance_indices.py EXP Y1 Y2

You can also run it from the specific python script in ``ecmean`` library ::

        ./performance_indices.py EXP Y1 Y2

- `Positional arguments`:

  EXP                   
    experiment identification

  Y1                    
    starting year of analysis

  Y2                   
    final year of analysis

- `Optional arguments`:

  -h, --help            show this help message and exit
  -s, --silent          do not print anything to std output
  -v LOGLEVEL, --loglevel LOGLEVEL
                        define the level of logging. default: error
  -j NUMPROC            number of processors to use
  -c CONFIG, --config CONFIG
                        set up a specific configuration file (config.yml is default)
  -i INTERFACE, --interface INTERFACE
                        set up a specific interface file (override config.yml)
  -m MODEL, --model MODEL
                        model name
  -e ENSEMBLE, --ensemble ENSEMBLE
                        variant label (i.e. ripf number for cmor)
  -d, --debug           activate cdo debugging
  -k CLIMATOLOGY        which climatology you want to use (EC23: default, RK08 alternative under update)
  -r RESOLUTION         only EC22: resolution of the climatology (r180x90 or r360x180)

Example 
^^^^^^^

Usage example for CMIP6 (running on 12 cores for EC-Earth3 historical)::

  > .performance_indices historical 1990 1999 -j 12 -m EC-Earth3 -e r1i1p1f1 -i CMIP6 

Usage example for EC-Earth4 (running on 4 cores for EC-Earth4 experment ABC1)::

  > performance_indices ABC1 1990 1999 -j 4


Output
^^^^^^

The result is produced in a form a YAML file, indicating PIs for each variable, region and season, that can be stored for later evaluation. 
Most importantly, a figure is produced showing a score card for the different regions, variables and seasons.
For the sake of simplicity, the PIs figure is computed as the ratio between the model PI and the average value estimated over the (precomputed) ensemble of CMIP6 models. 
An example of the the output for a single year of the EC-Earth3 historical simulation is shown here below.

.. figure:: _static/pitestfigure.png
   :align: center
   :width: 600px
   :alt: PI for ECearth3

   An example for a single year of the EC-Earth3 historical r1i1p1f1 simulation.

.. note ::
  Such plot is currently available for the EC23 climatology only, which is currently computed on a 30-year time window from 1990 to 2019 using about 10 models.
  Similarly, season-dependent computation are available only for EC23.
  Details on the field used are reported in the ``ecmean/climatology/EC23/pi_climatology_EC23.yml`` file.


Climatologies available
^^^^^^^^^^^^^^^^^^^^^^^

The performance indices built on the comparison between model data and a pre-computed climatology of several variables.
The ECmean climatology - from the previous CDO-based code - is currently defined as ``RK08``, and although still available, is not reccmmented for use since it is based on old datasets (e.g. ERA40). 

A new climatology has been developed making use of high-resolution data (e.g. CRU, ERA5, MSWEP, etc.) and is now defined as ``EC23``, using a 1x1 deg resolution and being the new default. 
An intermediate version knowns as ``EC22`` is available but not recommended and will be removed soon.

Properties of each climatology - as which variables, dataset and CMIP6 models used - can be inspected looking at ``ecmean/climatology/{clim}/pi_climatology_{clim}.yml`` files.

Climatology computation
^^^^^^^^^^^^^^^^^^^^^^^

Climatology is computed by the ``ecmean/utils/clim-create.py`` script, which is included in the repository for documentation.
It is based on a YML file which is tells the script where to retrieve the data, identifying all the required data folder and names. 
The tool loops over the variable and produces the yearly and seasonal average of the climate, as well as the interannual variance required for PIs. 
To avoid that grid points with unrealistic low variance affect the computation of the PIs, a filter to exclude outlier is introduced. This is based on the 5-sigma of the log10 distribution of each variable and each season.

Once the climatology is created, the script ``ecmean/utils/cmip6-clim-evaluate.py`` is used to run iteratively on a set of CMIP6 models and to compute the multi model mean of the PIs (for each region and season).
This is then stored in the ``ecmean/climatology/{clim}/pi_climatology_{clim}.yml`` and then used to provide a ratio between the original PI and the CMIP6 ensemble. 

