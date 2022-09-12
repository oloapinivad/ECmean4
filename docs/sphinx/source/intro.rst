Introduction
============

ECmean4 is a lightweight parallelized tool for evaluation of basic properties of Global Climate Models, as global mean, climate model performance indices and pattern correlation (under develpoment).
 
It builts on the original `ECmean <https://github.com/plesager/ece3-postproc/tree/master/ECmean>`_ which has been used for EC-Earth2 and EC-Earth3 evaluation, but it uses Python3 use YML configuration files. 
While the original ECmean4 version has been developed via CDO lazy calls, the current version is based on Xarray+Dask.


Under the hood
--------------

ECmean4 is built on a `Xarray <https://docs.xarray.dev/en/stable/>`_ + `Dask <https://examples.dask.org/xarray.html>`_ lazy calls which are executed in a single instance at the end of the script, exploiting paralellization on multiple variables. 
This allows to have a fast data analysis without writing unnecessary files on disk. Interpolation is carried out with `xESMF <https://xesmf.readthedocs.io/en/latest/>`_. 
Working with YML files allows for a more flexibile usage of two different climatologies and makes possible adding new climate models once a new interface file is developed. 

In order to assess that all the computations are correct, masks and area weights are pre-computed to provide robust global averages for all the different variables. 

ECmean4 also takes into account possible unit mismatch between the original dataset and the observational datasets, making using of the MetPY extension of the Pint python package. Heat and moisture flux sign convention is also assessed.

For the performance indices, since interpolation is required, weights are pre-computed only once to increase efficiency. Although conservative interpolation would be the better option, so far bilinear interpolation is preferred since it ensures more consistent results. 

.. note ::
	The original code has been developed using on CDO, but it has been replaced with xarray due to computational and scalability reasons.
	The new code is completely backward compatible, although some small differences in the computation are found due to interpolation. The older code is still available in the ``cdo`` subfolders. 
	
Computational performances
--------------------------

This has been tested on a single login node of ECMWF Atos AA making use of 10 years of a coupled EC-Earth4 run at Tco95L95-ORCA1

.. list-table:: ECmean4 single core performance
   :widths: 25 25 25
   :header-rows: 1

   * - SCRIPT
     - CDO version
     - Xarray version
   * - Global Mean
     - 1m39s
     - 37s
   * - Performance Indices
     - 2m2s
     - 32s

Scaling with multiple cores has been shown to be better for global mean than performance indices. 

.. note ::
	So far we cannot exploit completely dask in the xarray version due to previous existence of the multiprocessing library which seems to be conflicting. 
	The scheduler is thus set to synchronous with ``dask.config.set(scheduler="synchronous")``
