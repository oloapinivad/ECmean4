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
        self.ATMGRIDFILE = ''
        self.OCEGRIDFILE = ''
        self.LMFILE = ''
        self.SMFILE = ''
        self.ATMGAFILE = ''
        self.OCEGAFILE = ''
        self.ATMFIX = ''
        self.OCEFIX = ''
        self.TMPDIR = tempdir
        self.tempstore = CdoTempfileStore(dir=tempdir)
        self.cdo = Cdo(*args, **kwargs)
        self.domain = ''
        self.infile = ''

    def set_grids(self, atminifile, oceinifile): 
        """Load grid description files for both atmosphere and ocean"""
       
        self.ATMGRIDFILE = self.tempstore.newFile()
        griddes = self.cdo.griddes(input=str(atminifile))
        with open(self.ATMGRIDFILE, 'w', encoding='utf-8') as f:
            for line in griddes:
                print(line, file=f)

        self.OCEGRIDFILE = self.tempstore.newFile()
        griddes = self.cdo.griddes(input=str(oceinifile))
        with open(self.OCEGRIDFILE, 'w', encoding='utf-8') as f:
            for line in griddes:
                print(line, file=f)

    def make_masks(self, atminifile, extra = '') : 
        """Create land-sea masks for atmosphere model"""
    
        # prepare ATM LSM: this need to be improved, since it is clearly model dependent
        atmfix = eval(f"f'{self.ATMFIX}'")
        self.LMFILE = self.cdo.selname('LSM',
                                       input=f'-setctomiss,0 -gec,0.5 {extra} {atmfix} {atminifile}',
                                       options='-t ecmwf -f nc')
        self.SMFILE = self.cdo.addc('1', input=f'-setctomiss,1 -setmisstoc,0 {self.LMFILE}')

    def make_grids_areas(self, atminifile, oceinifile) : 
        """Create grid area file for initial files"""
    
        # this need to be improved, since it is clearly model dependent
        atmfix = eval(f"f'{self.ATMFIX}'")
        self.ATMGAFILE = self.cdo.gridarea(input=f'{atmfix} {atminifile}')
        self.OCEGAFILE = self.cdo.expr('area=e1t*e2t', input=oceinifile)

    def set_fixgrid(self, atmdomain, ocedomain) :
        """Define the command require for correcting model grid"""
        """String should be evaluated later, when corresponding files are defined"""
        
        # this could improved using the modelname variable: if EC-Earth, do this...
        if atmdomain == 'oifs' :
            self.ATMFIX = '-setgridtype,regular -setgrid,{self.ATMGRIDFILE}'
        else :
            sys.exit('Atmospheric component not supported')
        if ocedomain == 'nemo' :
            self.OCEFIX = '-setgridarea,{self.OCEGAFILE}'
        else :
            sys.exit('Oceanic component not supported')
            

    #def make_grids(self, atminifile, oceinifile, extra=''):
    #    """Initialize some useful helper files"""
    #
    #    # prepare ATM LSM
    #
    #    fix = f'{extra} -setgridtype,regular -setgrid,{self.ATMGRIDFILE}'
    #    self.LMFILE = self.cdo.selname('LSM',
    #                                   input=f'-setctomiss,0 -gec,0.5 {fix} {atminifile}',
    #                                   options='-t ecmwf -f nc')
    #    self.SMFILE = self.cdo.addc('1', input=f'-setctomiss,1 -setmisstoc,0 {self.LMFILE}')
    #    self.ATMGAFILE = self.cdo.gridarea(input=f'{self.LMFILE}')
    #
    #    # prepare OCE areas
    #    self.OCEGAFILE = self.cdo.expr('area=e1t*e2t', input=oceinifile)

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

        if not self.domain:
            sys.exit('Needed to define a domain with setdomain() method first')

        # this should be replaced for a more general "ocean" or "atmosphere"
        # currently using eval to implement string formatting, there should a better way
        if self.domain == 'nemo':
            ocefix = eval(f"f'{self.OCEFIX}'")
            self.pipe = ocefix +  ' ' + self.pipe
        elif self.domain == 'oifs':
            atmfix = eval(f"f'{self.ATMFIX}'")
            self.pipe = atmfix + ' ' + self.pipe

    def mask(self, mask_type):
        if not self.LMFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')

        if mask_type == 'land':
            self.mul(self.LMFILE)
        elif mask_type in ['sea', 'ocean']:
            self.mul(self.SMFILE)

    def masked_meansum(self, mask_type):
        if not self.LMFILE:
            sys.exit('Needed grid file not defined, call make_grids method first')

        if mask_type == 'land':
            self.chain(f'fldsum -mul {self.ATMGAFILE} -mul {self.LMFILE}')
        elif mask_type in ['sea', 'ocean']:
            self.chain(f'fldsum -mul {self.ATMGAFILE} -mul {self.SMFILE}')
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
            self.start()  # clear pipe
        return out

    def output(self, infile, **kwargs):
        return float(self.execute('output', input=infile, **kwargs)[0])
