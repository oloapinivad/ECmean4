Configuration
=============

Configuration file
------------------
A template configuration file is included in the repository, including the folder definition and all the details that can be manually adjusted. 
You will need to copy the original template to create your own local configuration file :

.. code-block:: shell

    > cp config.tmpl config.yml 

and then edit the required folder before running the script. 

interface
	The interface to access the file (see below). So far supported interfaces are `EC-Earth4`, `CMIP6` and `CMIP5`. 
model	
	The mode name: it could be `EC-Earth4` if running with the same interface, but could be any of the CMIP5/CMIP6 models.
dirs: exp
	Where the data from the experiments are
dirs: tab
	Where the output shoul be placed
dirs: clm
	Where the ECmean4 climatology is installed, i.e. the `ECmean4/climatology` folder

The configuration files defines also the variables which are assessed for the global fields an for the PIs, divided between oceanic and atmospheric variables. This should not be modifieds. 

Interface files
---------------

Conversion from model variables - as well as the correspondent file structure - is handled by a interface-depedent file ``interfaces/interface_*.yml`` that can support EC-Earth4, CMIP5 and CMIP6 models. 

.. note::
	It is not necessary to modify this, but it could be required if - for example - your CMIP5/6 directory tree does not reflect exactly the one available on ESGF. 

