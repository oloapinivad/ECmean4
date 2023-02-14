Configuration
=============

Configuration file
------------------
A template configuration file is included in the repository, including the folder definition and all the details that can be manually adjusted. 
You will need to copy the original template to create your own local configuration file ::
	
    > cp config.tmpl config.yml 

and then edit the required folder before running the script. 

General configuration
---------------------

These are the mandatory properties to be set up before running ECmean4.

interface
	The interface to access the model files (see below). So far supported interfaces are ``EC-Earth4``, ``CMIP6``, ``CMIP5``. 
model	
	The mode name: it could be ``EC-Earth4`` if running with the same interface, but could be any of the CMIP5/CMIP6 models.
dirs: exp
	Where the data from the experiments are
dirs: tab
	Where the output should be placed
dirs: fig
	Where the output figures should be placed (for ``performance_indices`` only)
dirs: clm
	Where the ECmean4 climatology is installed, i.e. the ``ECmean4/climatology`` folder

.. note::
	You could call a specific configuration file with ``-c config_myconfig.yml`` when executing the ECmean4 commands (see Usage section), so that you can have multiple models on the same machine.

Global Mean configuration
-------------------------

The configuration files defines also the details of the global mean. 
Modification to this is not necessary, but please be aware that the variable convention follow a CMOR-like definition.

atm_vars: 
	the list of atmospheric fields to be processed

oce_vars: 
	the list of oceanic fields to be processed

ice_vars:
        the list of sea-ice fields to be processed

tab_vars: 
	the list of fields to be reported in the output table

Performance indices configuration
---------------------------------

The configuration files defines also the details of the performance indices. 
Modification to this is not necessary, but please be aware that the variable convention follow a CMOR-like definition.

2d_vars: 
	the list of atmospheric 2d fields to be processed

3d_vars: 
	the list of atmospheric 3d fields to be processed. For this, PIs are computed based on the zonal mean cross-section.

oce_vars: 
	the list of oceanic 2d fields to be processed

ice_vars: 
	the list of sea-ice 2d fields to be processed

regions: 
	the list of regions on which compute the PI. Four regions are supported. Supported regions are: ``Global`` (90S-90N), ``North Midlat`` (30N-90N), ``Tropical`` (30S-30N), ``South Midlat`` (90S-30S) 

seasons:
	the list of seasons on which compute the PI. Four standard seasons are supported expressed as 3-string character (e.g. ``DJF``). ``ALL`` defines the yearly average. Default includes yearly, winter and summer.

resolution:
	the resolution on which compute PI. Do not change. 


Interface files
---------------

Conversion from model variables - as well as the correspondent file structure - is handled by a interface-dependent file ``ecmean/interfaces/interface_*.yml`` that can support EC-Earth4, CMIP5 and CMIP6 models. 
New interface files can be developed exploiting of the flexible file handling increasing the range of supported models. 

It is important to provide land-sea mask for both atmospheric and oceanic grid, and when possible - recommended for non-regular grid - it is important to provide also files indicating grid cell areas.
These latter are used for interpolation and weighted averages. ECmean4 has a few routines that tries to compute each grid cell area from lon/lat boundaries, but this can fail for non-regular grids.

.. note::
	It is not necessary to modify the interface file, but it could be required if - for example - your CMIP5/6 directory tree does not reflect exactly the one available on ESGF. 

Of course, if some specific management of the model grid is requested - for example if the model grid does not follow standard naming convention - it might be necessary to intervene directly within the code. 
A function named ``xr_preproc()``  within ``ecmean/libs/ncfixers.py``  might be adjusted in the case to allow for a correct interpretation of the model grid.


CMOR compatibility
------------------

It is possible to use ECmean4 tools also to analyze CMOR-like files for CMIP5 or CMIP6. This assumes a standard ESGF directory structure but you can change it by modifying the corresponding interface files ``ecmean/interfaces/interface_CMIP6.yml`` and ``ecmean/interfaces/interface_CMIP6.yml``.
In order to allow masking and interpolation you will need the ``sftlf`` (mandatory), ``sftof`` (not required, but suggested) and ``areacello`` (mandatory) variables for you experiment of interest too.


