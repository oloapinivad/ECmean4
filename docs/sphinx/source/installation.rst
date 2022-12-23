Installation
============

Getting the code
----------------

You can freely clone from the Github Repository ::

    > git clone https://github.com/oloapinivad/ECmean4.git
    
.. note ::

    Please note that if you clone with HTTPS you will not be able to contribute to the code, even if you are listed as collaborator.
    If you want to be a developer you should clone with SSH and you should add your own SSH key on the GitHub portal: 
    please check the `procedure on the Github website <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_ .


Install with Conda environment
------------------------------

It is recommended to work into a Conda environment with a recent Python3 version, which can be created with the proper environment file:
Once you set up Conda, ECmean4 can be easily installed with ::

    > conda env create --name ECmean4 -f environment.yml

with ``ECmean4`` is just an arbitrary name. Then you can activate the environment with ::

    > conda activate ECmean4


Requirements
------------

The required packages are listed in ``environment.yml``. 
A secondary environment available in  ``dev_environment.yml`` can be used for development. 

.. warning::
	Python >=3.8 is requested, but Python 3.11 is so far not supported due to conflicting packages




