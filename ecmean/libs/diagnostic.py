#!/usr/bin/env python3
'''
Shared diagnostic class for ECmean4
'''

import os
import logging
from pathlib import Path
import xarray as xr
import sys
from ecmean.libs.files import load_yaml

####################
# DIAGNOSTIC CLASS #
####################


class Diagnostic():
    """General container class for common variables"""

    def __init__(self, args):

        # arguments from command line/function
        self.expname = args.exp
        self.year1 = args.year1
        self.year2 = args.year2
        self.years_joined = list(range(self.year1, self.year2 + 1))
        self.fverb = not args.silent
        self.ftable = getattr(args, 'line', False)
        self.ftrend = getattr(args, 'trend', False)
        self.linefile = getattr(args, 'output', None)
        self.debug = getattr(args, 'debug', False)
        self.numproc = args.numproc
        self.modelname = getattr(args, 'model', '')
        self.climatology = getattr(args, 'climatology', 'EC23')
        self.interface = getattr(args, 'interface', '')
        self.resolution = getattr(args, 'resolution', '')
        self.ensemble = getattr(args, 'ensemble', 'r1i1p1f1')
        self.funcname = args.funcname.split(".")[1]
        if self.year1 == self.year2:
            self.ftrend = False
        print(self.funcname)

        #  These are here in prevision of future expansion to CMOR
        self.grid = '*'
        self.version = '*'
        self.frequency = '*mon'

        # current path
        self.indir = Path(os.path.dirname(os.path.abspath(__file__)))

        # get the config file
        if args.config:
            cfg = load_yaml(args.config)
        else:
            cfg = load_yaml(self.indir / '../../config.yml')

        # Various input and output directories
        self.ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']))
        self.TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
        self.FIGDIR = Path(os.path.expandvars(cfg['dirs']['fig']))

        # init for global mean
        if self.funcname == 'global_mean':
            self.cfg_global_mean(cfg)

        # init for performance indices
        if self.funcname in 'performance_indices':
            self.cfg_performance_indices(cfg)

        # setting up interface file
        if not self.interface:
            self.interface = cfg['interface']
        if not self.modelname:
            self.modelname = cfg['model']['name']

        # allow for both interface name or interface file
        if not os.path.exists(self.interface):
            self.interface = self.indir / Path(
                '../interfaces',
                f'interface_{self.interface}.yml')

        # load the possible xarray dataset
        if args.xdataset is not None:
            if isinstance(args.xdataset, (xr.DataArray, xr.Dataset)):
                logging.warning('You asked to use your own xarray dataset/datarray...')
                self.xdataset = args.xdataset
            else:
                sys.exit('Cannot used the xdataset, is not Xarray object')
        else:
            self.xdataset = None

    def cfg_global_mean(self, cfg):
        """Set up configuration details for global mean"""

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

        if self.linefile:
            self.ftable = True
        else:
            self.linefile = self.TABDIR / 'global_means.txt'

    def cfg_performance_indices(self, cfg):
        """Set up configuration for performance indices"""

        self.regions = cfg['PI']['regions']
        self.seasons = cfg['PI']['seasons']
        self.field_2d = cfg['PI']['2d_vars']['field']
        self.field_3d = cfg['PI']['3d_vars']['field']
        self.field_oce = cfg['PI']['oce_vars']['field']
        self.field_ice = cfg['PI']['ice_vars']['field']
        self.field_all = self.field_2d + self.field_3d + self.field_oce + self.field_ice

        # hard-coded resolution (due to climatological dataset)
        if self.climatology == 'RK08':
            logging.error('RK08 can work only with r180x91 grid')
            self.resolution = 'r180x91'
        else:
            if not self.resolution:
                self.resolution = cfg['PI']['resolution']

        # hard-coded seasons (due to climatological dataset)
        if self.climatology in ['EC22', 'RK08']:
            logging.error('only EC23 climatology support multiple seasons! Keeping only yearly seasons!')
            self.seasons = ['ALL']

        self.CLMDIR = Path(
            os.path.expandvars(
                cfg['dirs']['clm']),
            self.climatology)
        self.RESCLMDIR = Path(self.CLMDIR, self.resolution)
        self.climfile = self.CLMDIR / f'pi_climatology_{self.climatology}.yml'