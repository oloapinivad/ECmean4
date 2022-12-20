Configuration
=============

Configuration file
------------------
A template configuration file is included in the repository, including the folder definition and all the details that can be manually adjusted. 
You will need to copy the original template to create your own local configuration file :

.. code-block:: shell

    > cp config.tmpl config.yml 

and then edit the required folder before running the script. 

General configuration
---------------------

These are the mandatory properties to be set up before running Ecmean4.

interface
	The interface to access the file (see below). So far supported interfaces are ``EC-Earth4``, ``CMIP6``, ``CMIP5``. 
model	
	The mode name: it could be `EC-Earth4` if running with the same interface, but could be any of the CMIP5/CMIP6 models.
dirs: exp
	Where the data from the experiments are
dirs: tab
	Where the output shoul be placed
dirs: tab
	Where the output figures should be placed (for ``performance_indices.py`` only)
dirs: clm
	Where the ECmean4 climatology is installed, i.e. the ``ECmean4/climatology`` folder

The configuration files defines also the details of the global mean and of the performance indices. Modification to this is not necessary. 

.. note::
	You could call a specific configuration file with ``-c config_myconfig.yml`` when executing the ECmean4 command (see Usage section), so that you can have multiple models on the same machine.

Global Mean configuration
-------------------------

The configuration files defines also the details of the global mean. Modification to this is not necessary. 

atm_vars: 
	the list of atmospheric fields to be processed

oce_vars: 
	the list of oceanic fields to be processed

tab_vars: 
	the list of fields to be reported in the output table

Performance indices configuration
---------------------------------

The configuration files defines also the details of the performance indices. Modification to this is not necessary. 

2d_vars: 
	the list of atmospheric 2d fields to be processed

3d_vars: 
	the list of atmospheric 3d fields to be processed. For this, PI are computed based on the zonal mean cross-section.

oce_vars: 
	the list of oceanic 2d fields to be processed

ice_vars: 
	the list of sea-ice 2d fields to be processed

regions: 
	the list of regions on which compute the PI. Four regions are supported. Please see the table below.

seasons:
	the list of seasons on which compute the PI. Four standard seasons are supported expressed as 3-string character (e.g. DJF). ``ALL`` defines the yearly average.

resolution:
	the resolution on which compute PI. Do not change. 


Interface files
---------------

Conversion from model variables - as well as the correspondent file structure - is handled by a interface-depedent file ``interfaces/interface_*.yml`` that can support EC-Earth4, CMIP5 and CMIP6 models. 
New interface files can be developed exploiting of the flexible file handling increasing the range of supported models.

.. note::
	It is not necessary to modify this, but it could be required if - for example - your CMIP5/6 directory tree does not reflect exactly the one available on ESGF. 


CMOR compatibility
------------------

It is possible to use these tools also to analyze CMOR files for CMIP5 or CMIP6. This assumes a standard ESGF directory structure but you can change it by modifying the corresponding interface files ``interfaces/interface_CMIP6.yml`` and ``interfaces/interface_CMIP6.yml``.
In order to allow masking you will need the `sftlf`, `sftof` and `areacello` variables for you experiment of interest too.


