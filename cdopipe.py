#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 Class to create pipes of cdo commands

 @author Jost von Hardenberg (jost.hardenberg@polito.it), March 2022
'''

from cdo import Cdo, CdoTempfileStore
import tempfile

cdo = Cdo()

class CdoPipe:
    """A class to add commands in sequence to a cdo pipe"""

    def __init__(self, tempdir=tempfile.gettempdir()):
        """Initialize object pipe"""
        self.pipe = ''
        self.GRIDFILE = ''
        self.LMFILE = ''
        self.SMFILE = ''
        self.GAFILE = ''
        self.OCEGAFILE = ''
        self.TMPDIR = tempdir
        self.tempStore = CdoTempfileStore(dir = tempdir)

    def make_grids(self, atminifile, oceinifile):
        """Initialize some useful helper files"""

        self.GRIDFILE = self.tempStore.newFile()
        #self.GRIDFILE=str(self.TMPDIR / f'grid.txt')
        griddes = cdo.griddes(input=atminifile)
        with open(self.GRIDFILE, 'w') as f:
            for line in griddes:
                print(line, file=f)

        # prepare ATM LSM
        self.LMFILE = cdo.selname('LSM', input=f'-setgridtype,regular {atminifile}', options='-t ecmwf')
        self.SMFILE = cdo.mulc('-1', input=f'-subc,1 {self.LMFILE}')
        self.GAFILE = cdo.gridarea(input=f'-setgridtype,regular {self.LMFILE}')

        # prepare OCE areas
        self.OCEGAFILE = cdo.expr('area=e1t*e2t', input=oceinifile)

    def chain(self, cmd):
        """Adds a generic cdo operator"""
        self.pipe = '-' + cmd + ' ' + self.pipe

    def start(self, domain):
        """Cleans pipe for a new application, requires specifying the domain"""
        # ocean variables require specifying grid areas
        # atm variables require fixing the grid
        if not self.GRIDFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')

        if(domain=='oce'):
            self.pipe = f'-setgridarea,' + self.OCEGAFILE
        else:
            self.pipe = f'-setgridtype,regular -setgrid,{self.GRIDFILE}'

    def masked_mean(self, mask_type):
        if not self.LMFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')

        if mask_type == 'land':
            self.chain(f'fldsum -mul {self.GAFILE} -mul {self.LMFILE}')
        elif mask_type in ['sea', 'ocean']:
            self.chain(f'fldsum -mul {self.GAFILE} -mul {self.SMFILE}')
        else:
            self.chain('fldmean')

    def selname(self, var):
        self.chain(f'selname,{var}')

    def expr(self, var, expr):
        self.chain(f'expr,{var}={expr}')

    def selcode(self, var):
        self.chain(f'selcode,{var}')

    def timmean(self):
        self.chain(f'timmean')

    def output(self, infile):
        cmd = self.pipe + f' {infile}'
        return float(cdo.output(input=cmd)[0])

    def execute(self, *args, **kwargs):
        return cdo(self.pipe, *args, **kwargs)

