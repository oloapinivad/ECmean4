Installation
============

Prerequisites
-------------

Python 3.7+ is required. 

Getting the code
----------------

You can freely clone from the Github Repository ::

    > git clone git@github.com:oloapinivad/ECmean4.git


Xarray version: Install with Conda environment
----------------------------------------------

It is recommended to work into a Conda environment with a recent Python3 version, which can be created with the proper environment file:
Once you set up Conda, ECmean4 can be easily installed with:

.. code-block:: shell

    > conda env --name ECmean4 create -f docs/environment.yml

with ``ECmean4`` is just an arbitrary name. Then you can activate the environment with:

.. code-block:: shell

    > conda activate ECmean4

Requirements
------------

The required packages are listed in `docs/environment.yml`` 


CDO version: Install in a python virtual environment
----------------------------------------------------

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



