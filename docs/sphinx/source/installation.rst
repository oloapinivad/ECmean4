Installation
============

Prerequisites
-------------

Python 3.7+ is required. 

Getting the code
----------------

You can freely clone from the Github Repository ::

    > git clone git@github.com:oloapinivad/ECmean4.git

.. note ::

    To be able to use the SSH access to GitHub, you shoud add your own ssh key: 
    please check the `procedure on the Github website <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_
    Alternatively you can clone with HTTPS, but you will not be able to push anything to the repo.
    


Xarray version: Install with Conda environment
----------------------------------------------

It is recommended to work into a Conda environment with a recent Python3 version, which can be created with the proper environment file:
Once you set up Conda, ECmean4 can be easily installed with:

.. code-block:: shell

    > conda env --name ECmean4 create -f environment.yml

with ``ECmean4`` is just an arbitrary name. Then you can activate the environment with:

.. code-block:: shell

    > conda activate ECmean4

Requirements
------------

The required packages are listed in `environment.yml`` 



