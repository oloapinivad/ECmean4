#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import os
import logging
from pathlib import Path
import pandas as pd
import re


####################
# DIAGNOSTIC CLASS #
####################


class Diagnostic():
    """General container class for common variables"""

    def __init__(self, args, cfg):
        self.expname = args.exp
        self.year1 = args.year1
        self.year2 = args.year2
        self.fverb = not args.silent
        self.ftable = getattr(args, 'line', False)
        self.ftrend = getattr(args, 'trend', False)
        self.debug = getattr(args, 'debug', False)
        self.numproc = args.numproc
        self.modelname = getattr(args, 'model', '')
        self.climatology = getattr(args, 'climatology', 'EC23')
        self.interface = getattr(args, 'interface', '')
        self.resolution = getattr(args, 'resolution', '')
        self.regions = cfg['PI']['regions']
        self.seasons = cfg['PI']['seasons']
        if not self.modelname:
            self.modelname = cfg['model']['name']
        if self.year1 == self.year2:  # Ignore if only one year requested
            self.ftrend = False
        #  These are here in prevision of future expansion to CMOR
        if not self.interface:
            self.interface = cfg['interface']
        self.frequency = '*mon'
        self.ensemble = getattr(args, 'ensemble', 'r1i1p1f1')
        self.grid = '*'
        self.version = '*'

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

        # Various input and output directories
        self.ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']))
        self.TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
        self.FIGDIR = Path(os.path.expandvars(cfg['dirs']['fig']))
        self.CLMDIR = Path(
            os.path.expandvars(
                cfg['dirs']['clm']),
            self.climatology)
        self.RESCLMDIR = Path(self.CLMDIR, self.resolution)
        self.years_joined = list(range(self.year1, self.year2 + 1))

        if hasattr(args, 'output') and args.output:
            self.linefile = args.output
            self.ftable = True
        else:
            self.linefile = self.TABDIR / 'global_means.txt'


##################
# HELP FUNCTIONS #
##################


def get_variables_to_load(var, face):
    """Function to extract from the interface file the list of derived variable,
    i.e. the real variables to be loaded, for each of the cmorname introduced in the
    interface file

    Args:
        var: the cmorname variable of the data to be loaded
        face: the interface file
    """

    if 'derived' in face['variables'][var].keys():
        cmd = face['variables'][var]['derived']
        dervars = re.findall("[a-zA-Z0-9]+", cmd)
    else:
        dervars = [var]

    return dervars


def is_number(s):
    """Check if input is a float type"""

    try:
        float(s)
        return True
    except ValueError:
        return False


def numeric_loglevel(loglevel):
    """Define the logging level """
    # log level with logging
    # currently basic definition trought the text
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    return numeric_level


# def chunks(iterable, num):
#     """Generate num adjacent chunks of data from a list iterable
#     Split lists in a convenient way for a parallel process"""

#     size = int(np.ceil(len(iterable) / num))
#     it = iter(iterable)
#     return iter(lambda: tuple(itertools.islice(it, size)), ())

# def split(iterable, num):
#     k, m = divmod(len(iterable), num)
#     return (iterable[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(num))

def runtime_weights(varlist):
    """Define the weights to estimate the best repartition of the cores
    This is done a-priori, considering that 1) compound variables are more difficult to compute
    2) 3d variables requires more evaluation"""

    w = {}
    for k in varlist:
        if k in ['ua', 'ta', 'va', 'hus']:
            t = 8
        elif k in ['pme', 'net_sfc_nosn', 'net_sfc', 'toamsfc_nosn', 'toamsfc',
                   'pr_oce', 'pme_oce', 'pr_land', 'pme_land', 'net_toa']:
            t = 3
        else:
            t = 1
        w[k] = t

    return w


def weight_split(a, n):
    """use the weights by runtime_weights to provide a number of n chunks based on the
    minimum possible value of each chunk computing time. Then, provide the list of variables
    to be used in the multiprocessing routine"""

    ws = sorted(runtime_weights(a).items(), key=lambda x: x[1], reverse=True)
    ordered = dict(ws)

    elists = [[0] for _ in range(n)]
    olists = [[] for _ in range(n)]

    count = 0
    for f in ordered.keys():
        elists[count].append(ordered[f])
        olists[count].append(f)
        count = elists.index(min(elists))

    return olists


def get_domain(var, face):
    """Given a variable var extract its domain (oce or atm) from the interface.
    To do so it creates a dictionary providing the domain associated with a component.
    (the interface file specifies the component for each domain instead)"""

    comp = face['filetype'][face['variables'][var]['filetype']]['component']
    d = face['model']['component']
    domain = dict(zip([list(d.values())[x]
                  for x in range(len(d.values()))], d.keys()))
    return domain[comp]


# def get_component(face):  # unused function
#     """Return a dictionary providing the domain associated with a variable
#     (the interface file specifies the domain for each component instead)"""

#     d = face['component']
#     p = dict(zip([list(d.values())[x]['domain']
#              for x in range(len(d.values()))], d.keys()))
#     return p


####################
# OUTPUT FUNCTIONS #
####################

def dict_to_dataframe(varstat):
    """very clumsy conversion of the nested 3-level dictionary
    to a pd.dataframe: NEED TO BE IMPROVED"""
    data_table = {}
    for i in varstat.keys():
        pippo = {}
        for outerKey, innerDict in varstat[i].items():
            for innerKey, values in innerDict.items():
                pippo[(outerKey, innerKey)] = values
        data_table[i] = pippo
    data_table = pd.DataFrame(data_table).T
    return data_table


def write_tuning_table(linefile, varmean, var_table, diag, ref):
    """Write results appending one line to a text file.
       Write a tuning table: need to fix reference to face/ref"""

    if not os.path.isfile(linefile):
        with open(linefile, 'w', encoding='utf-8') as f:
            print('%model  ens  exp from   to ', end='', file=f)
            for var in var_table:
                print('{:>12s}'.format(var), end=' ', file=f)
            print('\n%                         ', end=' ', file=f)
            for var in var_table:
                print('{:>12s}'.format(ref[var]['units']), end=' ', file=f)
            print(file=f)

    with open(linefile, 'a', encoding='utf-8') as f:
        print(f'{diag.modelname} {diag.ensemble} {diag.expname}',
              '{:4d} {:4d} '.format(diag.year1, diag.year2), end='', file=f)
        for var in var_table:
            print(
                '{:12.5f}'.format(
                    varmean[var] *
                    ref[var].get(
                        'factor',
                        1)),
                end=' ',
                file=f)
        print(file=f)
