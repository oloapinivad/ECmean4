API Reference
=============

This section provides a detailed reference for |ecmean|'s Application Programming Interface (API).

Usage from Python
-----------------

Both diagnostics are available as a function, which runs the full workflow (computation, YAML output and figure),
and as a class, which allows to run the steps separately and to access the results.

Global mean::

    from ecmean import GlobalMean, global_mean

    # full workflow
    global_mean("historical", 1990, 1999, config="config.yml", numproc=4, interface="CMIP6", model="EC-Earth3")

    # step by step
    gm = GlobalMean("historical", 1990, 1999, config="config.yml", numproc=4)
    gm.prepare()
    gm.run()
    gm.store()
    gm.plot(returnfig=True)
    gm.varmean  # dictionary with the results

Performance indices::

    from ecmean import PerformanceIndices, performance_indices

    # full workflow, the figure can be skipped
    performance_indices("historical", 1990, 1999, config="config.yml", numproc=4, plot=False)

    # step by step
    pi = PerformanceIndices("historical", 1990, 1999, config="config.yml", numproc=4)
    pi.prepare()
    pi.run()
    pi.store()
    pi.plot()
    pi.varstat  # dictionary with the results

Other arguments common to the classes and the functions are ``loglevel``, ``ensemble``, ``outputdir`` and ``title``.
``config`` can be either the path to a configuration file or a dictionary with the same structure.
``xdataset`` allows to pass an already opened ``xarray`` object instead of reading the files defined by the interface.
``GlobalMean`` also accepts ``reference``, ``trend``, ``line`` and ``addnan``, while ``PerformanceIndices`` accepts ``climatology``.

Reference
---------

.. automodule:: ecmean
    :members:
    :undoc-members:
    :show-inheritance:
