Introduction
============

ECmean4 is a lightweight tool for evaluation of basic properties of Global Climate Models. 
It builts on the original ECmean which has been used for EC-Earth2 and EC-Earth3 evaluation, but it uses Python3 as a scripting language to perform lazy calls to CDO.

Scripts
============
Two scripts takes care of all the computation and produces a txt table:

- ``global_mean.py``: computes the global averages for many dynamical and physical fields. It compares the output against a set of climatological values defined in ``gm_reference.yml``
- ``performance_indices.py``: computes Reichler and Kim Performance Indices. It work on a regular 2x2 grid and it is based on old climatological assessment present in the original ECmean. The climatology is defined by ``pi_climatology.yml`` script. Climatology for PI is still VERY outdated.

Configuration
=============
A template configuration file is included in the repo, including folder definition and all the details that can be manually adjusted You will need to copy ``cp config.tmpl config.yml`` and then edit the required folder before running the script
Conversion from model variables - as well as file structure - is handled by a model depedent file ``interrfaces/interface_*.yml`` that can support so far EC-Earth4, CMIP5 and CMIP6 models.



