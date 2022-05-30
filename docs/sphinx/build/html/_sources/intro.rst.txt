Introduction
============

ECmean4 is a lightweight parallelized tool for evaluation of basic properties of Global Climate Models, as global mean, pattern correlation and climate model performance indices.
 
It builts on the original `ECmean <https://github.com/plesager/ece3-postproc/tree/master/ECmean>`_ which has been used for EC-Earth2 and EC-Earth3 evaluation, but it uses Python3 as a scripting language to perform lazy calls to CDO and makes use of YML configuration files.


Under the hood
--------------

ECmean4 is built on a specific `cdopipe.py` class which is used to chain CDO operators serially as required and the execute them in a single or max two instances at the end of the script, exploiting paralellization on multiple variables. This allows to have a fast call to CDO without writing unnecessary files on disk and without opening NetCDF files in Python. Working with YML files allows for a more flexibile usage of different climatologies - currently under develpoment - and makes possible adding new climate models once a new interface file is developed. 

In order to assess that all the computation are correct, masks and area weights are pre-computed to provide robust global averages for all the different variables. 

ECmean4 also takes into account possible unit mismatch between the original dataset and the observational datasets, making using of the MetPY extension of the Pint python package. Heat and moisture flux sign convention is also assessed.

For the performance indices, since interpolation is required, weights are pre-computed only once to increase efficiency. Although conservative interpolation would be the better option, so far bilinear interpolation is preferred since it ensures more consistent results. 


Computational performances
--------------------------

.. note ::
	This section is draft


