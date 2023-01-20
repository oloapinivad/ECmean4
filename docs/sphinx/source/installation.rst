Installation
============


ECmean4 is a lightweight python package. We recommended to do it within a conda or better a mamba environment. 
For the moment 


Installation
------------

You can freely clone from the Github Repository ::

    > git clone https://github.com/oloapinivad/ECmean4.git
    
.. note ::

    Please note that if you clone with HTTPS you will not be able to contribute to the code, even if you are listed as collaborator.
    If you want to be a developer you should clone with SSH and you should add your own SSH key on the GitHub portal: 
    please check the `procedure on the Github website <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_ .

Then you can through the ECmean4 folder

    > cd ECmean4

and then you can set up the conda/mamba environment

    > conda env create --name ECmean4 -f environment.yml

Then you should activate the environment

    > conda activate ECmean4

Finally you can install the development package of ECmean4

    > pip install -e .

Checking everything is ok
-------------------------

From now on the two command line function of ECmean4 (``global_mean`` and ``performance_indices``) should be available in your environment
You can test by running in command line and you should get something like this: 

    > global_mean
    > usage: global_mean [-h] [-s] [-t] [-l] [-o FILE] [-m MODEL] [-c CONFIG] [-v LOGLEVEL] [-j NUMPROC] [-e ENSEMBLE] [-i INTERFACE] EXP Y1 Y2 
    > global_mean: error: the following arguments are required: EXP, Y1, Y2

You can also run tests by simply calling ``pytest`` from the ECmean4 folder

    > pytest

Requirements
------------

The required packages are listed in ``environment.yml``. 
A secondary environment available in  ``ecmean/utils/dev_environment.yml`` can be used for development. 

.. warning::
	Python >=3.8 is requested, but Python 3.11 is so far not supported due to conflicting packages




