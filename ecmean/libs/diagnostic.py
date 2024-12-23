#!/usr/bin/env python3
'''
Shared diagnostic class for ECmean4
'''

import os
import logging
from pathlib import Path
import xarray as xr
from ecmean.libs.files import load_yaml
from ecmean import __version__ as version

####################
# DIAGNOSTIC CLASS #
####################

loggy = logging.getLogger(__name__)


class Diagnostic():
    """General container class for common variables"""

    def __init__(self, exp, year1, year2, config, funcname,
                 line=False, trend=False,
                 resolution=None, ensemble='r1i1p1', addnan=False,
                 interface=None, modelname=None, outputdir=None,
                 xdataset=None, silent=None,
                 numproc=1, climatology='EC23'):
        """
        Initialize the Diagnostic instance.

        Args:
            args: Arguments from command line/function.
        """

        # arguments from command line/function
        self.expname = exp
        self.year1 = year1
        self.year2 = year2
        self.years_joined = list(range(self.year1, self.year2 + 1))
        self.ftable = line
        self.ftrend = trend
        self.numproc = numproc
        self.climatology = climatology
        self.silent = silent
        self.resolution = resolution
        self.ensemble = ensemble
        self.interface = interface
        self.funcname = funcname
        self.version = version
        self.addnan = addnan
        if self.year1 == self.year2:
            self.ftrend = False
        print(f'Welcome to ECmean4 v{self.version}: Running {self.funcname} with {self.numproc} cores!')

        #  These are here in prevision of future expansion to CMOR
        self.grid = '*'
        self.version = '*'
        self.frequency = '*mon'

        # current path
        self.indir = Path(os.path.dirname(os.path.abspath(__file__)))

        # get the config file, can be a string or a dictionary
        if config:
            if isinstance(config, dict):
                cfg = config
            elif isinstance(config, str):
                cfg = load_yaml(config)
            else:
                raise ValueError('Cannot load the config file')
        else:
            cfg = load_yaml(self.indir / '../../config.yml')

        # Various raise and input and output directories
        if not cfg['dirs']['exp']:
            raise ValueError('No experiment directory defined in config file')
        self.ecedir = Path(os.path.expandvars(cfg['dirs']['exp']))
        if outputdir is None:
            if not cfg['dirs']['tab']:
                raise ValueError('No table directory defined in config file')
            self.tabdir = Path(os.path.expandvars(cfg['dirs']['tab']))
            if not cfg['dirs']['fig']:
                raise ValueError('No figure directory defined in config file')
            self.figdir = Path(os.path.expandvars(cfg['dirs']['fig']))
        else:
            self.tabdir = Path(os.path.join(outputdir, 'YAML'))
            self.figdir = Path(os.path.join(outputdir, 'PDF'))

        # init for global mean
        if self.funcname == 'GlobalMean':
            self.cfg_global_mean(cfg)

        # init for performance indices
        if self.funcname in 'PerformanceIndices':
            self.cfg_performance_indices(cfg)

        # setting up interface file
        self.interface = interface or cfg['interface']

        # setting up model name
        self.modelname = modelname or cfg['model']['name']

        # allow for both interface name or interface file
        if not os.path.exists(self.interface):
            self.interface = self.indir / Path(
                '../interfaces',
                f'interface_{self.interface}.yml')

        # load the possible xarray dataset
        if xdataset is not None:
            if isinstance(xdataset, (xr.DataArray, xr.Dataset)):
                loggy.warning('You asked to use your own xarray dataset/datarray...')
                self.xdataset = xdataset
            else:
                raise ValueError('Cannot used the xdataset, is not Xarray object')
        else:
            self.xdataset = None

    def filenames(self, kind):
        """
        Return the filename for the output.
        """

        if self.funcname == 'GlobalMean':
            head = 'global_mean'
        elif self.funcname == 'PerformanceIndices':
            head = f'PI4_{self.climatology}'
        else:
            raise ValueError('Unknown function name')

        figurename = f'{head}_{self.expname}_{self.modelname}_{self.ensemble}_{self.year1}_{self.year2}.{kind}'

        if kind in ['yml', 'txt']:
            return self.tabdir / figurename
        if kind in ['pdf', 'png']:
            return self.figdir / figurename
        raise ValueError('Unknown file type')

    def cfg_global_mean(self, cfg):
        """
        Set up configuration details for global mean.

        Args:
            cfg: Configuration details.

        Returns:
            None
        """

        self.regions = cfg['global']['regions']
        self.seasons = cfg['global']['seasons']

        # define the different list required
        self.var_atm = cfg['global']['atm_vars']
        self.var_oce = cfg['global']['oce_vars']
        self.var_ice = cfg['global']['ice_vars']
        self.var_table = cfg['global']['tab_vars']
        self.var_all = list(dict.fromkeys(
            self.var_atm +
            self.var_table +
            self.var_oce +
            self.var_ice))

        self.reffile = self.indir / '../reference/gm_reference_EC23.yml'

        if self.ftable:
            self.linefile = self.tabdir / 'global_means.txt'

    def cfg_performance_indices(self, cfg):
        """
        Set up configuration for performance indices.

        Args:
            cfg: Configuration details.

        Returns:
            None
        """

        self.regions = cfg['PI']['regions']
        self.seasons = cfg['PI']['seasons']
        self.field_2d = cfg['PI']['2d_vars']['field']
        self.field_3d = cfg['PI']['3d_vars']['field']
        self.field_oce = cfg['PI']['oce_vars']['field']
        self.field_ice = cfg['PI']['ice_vars']['field']
        self.field_all = self.field_2d + self.field_3d + self.field_oce + self.field_ice

        # hard-coded resolution (due to climatological dataset)
        if self.climatology == 'RK08':
            loggy.error('RK08 can work only with r180x91 grid')
            self.resolution = 'r180x91'
        else:
            if not self.resolution:
                self.resolution = cfg['PI']['resolution']

        # hard-coded seasons (due to climatological dataset)
        if self.climatology in ['EC22', 'RK08']:
            loggy.error('Only EC23 climatology supports multiple seasons! Keeping only yearly seasons!')
            self.seasons = ['ALL']

        self.clmdir = Path(self.indir, '../climatology', self.climatology)
        self.resclmdir = Path(self.clmdir, self.resolution)
        self.climfile = self.clmdir / f'pi_climatology_{self.climatology}.yml'
