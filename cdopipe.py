#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 Class to create pipes of cdo commands

 @author Jost von Hardenberg (jost.hardenberg@polito.it), March 2022
 @author Paolo Davini (p.davini@isac.cnr.it), March 2022
'''

import sys
import tempfile
from cdo import Cdo, CdoTempfileStore

class CdoPipe:
    """A class to add commands in sequence to a cdo pipe"""

    def __init__(self, *args, tempdir=tempfile.gettempdir(), **kwargs):
        """Initialize object pipe"""
        self.pipe = ''
        self.GRIDFILE = ''
        self.LMFILE = ''
        self.SMFILE = ''
        self.GAFILE = ''
        self.OCEGAFILE = ''
        self.TMPDIR = tempdir
        self.tempStore = CdoTempfileStore(dir = tempdir)
        self.cdo = Cdo(*args, **kwargs)
        self.domain = ''
        self.infile = ''

    def make_grids(self, atminifile, oceinifile):
        """Initialize some useful helper files"""

        self.GRIDFILE = self.tempStore.newFile()
        #self.GRIDFILE=str(self.TMPDIR / f'grid.txt')

        griddes = self.cdo.griddes(input=str(atminifile))

        with open(self.GRIDFILE, 'w') as f:
            for line in griddes:
                print(line, file=f)

        # prepare ATM LSM
        fix = '-setgridtype,regular -setgrid,{self.GRIDFILE}'
        self.LMFILE = self.cdo.selname('LSM',
                                       input=f'-gec,0.5 {fix} {atminifile}',
                                       options='-t ecmwf')
        self.SMFILE = self.cdo.mulc('-1', input=f'-subc,1 {self.LMFILE}')
        self.GAFILE = self.cdo.gridarea(input=f'{self.LMFILE}')

        # prepare OCE areas
        self.OCEGAFILE = self.cdo.expr('area=e1t*e2t', input=oceinifile)

    def chain(self, cmd):
        """Adds a generic cdo operator"""
        self.pipe = '-' + cmd + ' ' + self.pipe

    def setdomain(self, domain):
        """Specify variable domain: used to set needed grid manipulations"""
        self.domain = domain

    def start(self):
        """Cleans pipe for a new application"""
        # ocean variables require specifying grid areas
        # atm variables require fixing the grid
        self.pipe = '{infile}'

    def fixgrid(self):
        """Applies grid fixes, requires specifying the domain"""
        # ocean variables require specifying grid areas
        # atm variables require fixing the grid

        if not self.GRIDFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')
        if not self.domain : 
            sys.exit('Needed to define a domain with setdomain() method first') 

        if self.domain=='nemo':
            self.pipe = f'-setgridarea,{self.OCEGAFILE} ' + self.pipe
        elif self.domain=='oifs':
            self.pipe = f'-setgridtype,regular -setgrid,{self.GRIDFILE} ' + self.pipe

    def masked_mean(self, mask_type):
        if not self.LMFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')

        if mask_type == 'land':
            self.chain(f'fldsum -mul {self.GAFILE} -mul {self.LMFILE}')
        elif mask_type in ['sea', 'ocean']:
            self.chain(f'fldsum -mul {self.GAFILE} -mul {self.SMFILE}')
        else:
            self.chain('fldmean')

    def setname(self, var):
        self.chain(f'setname,{var}')

    def selname(self, var):
        self.chain(f'selname,{var}')

    def selectname(self, var):
        self.chain(f'select,name={var}')

    def expr(self, var, expr):
        self.chain(f'expr,{var}={expr}')

    def selcode(self, var):
        self.chain(f'selcode,{var}')

    def timmean(self):
        self.chain('timmean')

    def invertlat(self):
        self.chain('invertlat')

    def zonmean(self):
        self.chain('zonmean')

    def sqr(self):
        self.chain('sqr')

    def sub(self, fname):
        self.pipe = '-sub ' + self.pipe + ' ' + fname

    def mul(self, fname):
        self.pipe = '-mul ' + self.pipe + ' ' + fname

    def div(self, fname):
        self.pipe = '-div ' + self.pipe + ' ' + fname

    def levels(self, infile):
        return list(map(float, self.cdo.showlevel(input=infile)[0].split()))

    def set_infile(self, infile):
        self.infile = infile

    def execute(self, cmd, *args, input='', keep=False, **kwargs):
        fn = getattr(self.cdo, cmd)
        if not input:
            input = self.infile
        # print("EXE ",self.pipe)
        out = fn(input=self.pipe.format(infile=input), *args, **kwargs)
        if not keep:
            self.start() # clear pipe
        return out

    def output(self, infile, **kwargs):
        return float(self.execute('output', input=infile, **kwargs)[0])
