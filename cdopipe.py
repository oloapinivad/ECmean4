#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 Class to create pipes of cdo commands

 @author Jost von Hardenberg (jost.hardenberg@polito.it), March 2022
'''

from cdo import Cdo, CdoTempfileStore
import tempfile
import os

class CdoPipe:
    """A class to add commands in sequence to a cdo pipe"""

    def __init__(self, tempdir=tempfile.gettempdir(), *args, **kwargs):
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

    def make_grids(self, atminifile, oceinifile):
        """Initialize some useful helper files"""

        self.GRIDFILE = self.tempStore.newFile()
        #self.GRIDFILE=str(self.TMPDIR / f'grid.txt')

        griddes = self.cdo.griddes(input=str(atminifile))

        with open(self.GRIDFILE, 'w') as f:
            for line in griddes:
                print(line, file=f)

        # prepare ATM LSM
        self.LMFILE = self.cdo.selname('LSM', input=f'-setgridtype,regular {atminifile}', options='-t ecmwf')
        self.SMFILE = self.cdo.mulc('-1', input=f'-subc,1 {self.LMFILE}')
        self.GAFILE = self.cdo.gridarea(input=f'-setgridtype,regular {self.LMFILE}')

        # prepare OCE areas
        self.OCEGAFILE = self.cdo.expr('area=e1t*e2t', input=oceinifile)

    def chain(self, cmd):
        """Adds a generic cdo operator"""
        self.pipe = '-' + cmd + ' ' + self.pipe

    def domain(self, domain):
        """Specify variable domain: used to set needed grid manipulations"""
        self.domain = domain

    def start(self, domain='', **kwargs):
        """Cleans pipe for a new application"""
        # ocean variables require specifying grid areas
        # atm variables require fixing the grid
        if domain:
            self.domain = domain
        self.pipe = '{infile}'

    def fixgrid(self, domain='', **kwargs):
        """Applies grid fixes, requires specifying the domain"""
        # ocean variables require specifying grid areas
        # atm variables require fixing the grid

        if not self.GRIDFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')

        if(self.domain=='oce'):
            self.pipe = f'-setgridarea,{self.OCEGAFILE} ' + self.pipe
        else:
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

#    def remapbil(self, res, infile):
#        cmd = self.pipe.format(infile=infile)
#        out = self.cdo.remapbil(res, input=cmd, output='/home/ccjh/ece4/tools/ECmean4/test.nc', debug=True)
#
#        print("Output of rmpb is:" , out)
#        print("it does exist?: ", os.path.isfile(out))
#        return out

    def sub(self, fname):
        self.pipe = '-sub ' + self.pipe + ' ' + fname

    def mul(self, fname):
        self.pipe = '-mul ' + self.pipe + ' ' + fname

    def div(self, fname):
        self.pipe = '-div ' + self.pipe + ' ' + fname

    def levels(self, infile):
        return list(map(float, self.cdo.showlevel(input=infile)[0].split()))

    def set_input(self, infile):
        self.infile = infile

    def execute(self, cmd, *args, input='', keep=False, **kwargs):
        fn = getattr(self.cdo, cmd) 
        if not input: input = self.infile
        # print("EXE ",self.pipe)
        out = fn(input=self.pipe.format(infile=input), *args, **kwargs)
        if not keep:
            self.start() # clear pipe
        return out

    def output(self, infile, **kwargs):
        return float(self.execute('output', input=infile, **kwargs)[0])

