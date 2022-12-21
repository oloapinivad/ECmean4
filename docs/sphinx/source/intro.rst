Introduction
============

About
-----

ECmean4 is a lightweight parallelized tool for evaluation of basic properties of Global Climate Models, as global mean and climate model performance indices.
 
It builts on the original `ECmean <https://github.com/plesager/ece3-postproc/tree/master/ECmean>`_ which has been used for EC-Earth2 and EC-Earth3 evaluation, but it uses Python3 use YML configuration files. 
While the original ECmean4 version has been developed via CDO lazy calls, the current version is based on Xarray+Dask.


Under the hood
--------------

ECmean4 is built on a `Xarray <https://docs.xarray.dev/en/stable/>`_ + `Dask <https://examples.dask.org/xarray.html>`_ lazy calls which are executed in a single instance at the end of the script, 
exploiting paralellization on multiple variables with `Multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_. 
This allows to have a fast data analysis without writing unnecessary files on disk. Interpolation is carried out with `xESMF <https://xesmf.readthedocs.io/en/latest/>`_. 
Working with YML files in each configuration aspects allows for a more flexibile usage, making possible expanding the support to new climate models or to include new reference climatologies. 

ECmean4 also takes into account possible unit mismatch between the original dataset and the observational datasets, making using of the `MetPY <https://unidata.github.io/MetPy/latest/index.html>`_ 
extension of the Pint python package. Heat and moisture flux sign convention is also assessed.

For the performance indices, since interpolation is required, weights are pre-computed only once to increase efficiency. 
Although conservative interpolation would be the better option, so far bilinear interpolation is preferred since it ensures more consistent results. 

.. note ::
	The original code has been developed using on CDO, but it has been replaced with xarray due to computational and scalability reasons.
	The new code is completely backward compatible, although some small differences in the computation are found due to interpolation. The older code is still available in the ``cdo`` subfolders. 
	
Computational performances
--------------------------

ECmean4 can process many years and multiple variables in less than 10 minutes (assuming that output is provided from monthly means). 
Performance indices are implicitly slower than global mean, but with a few cores available both can be completed in a couple of minutes.
Since parallelization is done along variables, it does not make sense (especially for performance indices) to use more than 6 cores. 

This has been tested on a single core machine for 30 years of a coupled EC-Earth3 CMIP6 historical run (i.e. TL255L91, about 0.7x0.7 deg), using the default ``config.yml`` (for performance indices evaluating on 3 seasons and 4 regions).

.. figure:: _static/benchmark.png
   :align: center
   :width: 600px
   :alt: Benchmark for EC-Earth3

   A multi-core benchmarking for Global Mean and Performance Indices over 10 year of CMIP6 EC-Earth3 data

.. note ::
	So far we cannot exploit completely dask in the xarray version due to previous existence of the multiprocessing library which is partially conflicting. This issue will be addressed in future release, but so far dask scheduler is thus set to synchronous with ``dask.config.set(scheduler="synchronous")``.
