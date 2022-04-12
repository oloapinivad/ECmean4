#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 Class to create pipes of cdo commands

 @author Jost von Hardenberg (jost.hardenberg@polito.it), March 2022
 @author Paolo Davini (p.davini@isac.cnr.it), March 2022
'''

import sys
import tempfile
import logging
from cdo import Cdo, CdoTempfileStore


class CdoPipe:
    """A class to add commands in sequence to a cdo pipe"""

    def __init__(self, *args, tempdir=tempfile.gettempdir(), **kwargs):
        """Initialize object pipe"""
        self.pipe = ''
        self.ATMGRID = ''
        self.OCEGRID = ''
        self.LANDMASK = ''
        self.SEAMASK = ''
        self.ATMGRIDAREA = 'dummy'
        self.OCEGRIDAREA = 'dummy'
        self.ATMWEIGHTS = ''
        self.OCEWEIGHTS = ''
        self.atmfix = ''
        self.ocefix = ''
        self.TMPDIR = tempdir
        self.tempstore = CdoTempfileStore(dir=tempdir)
        self.cdo = Cdo(*args, **kwargs)
        self.domain = ''
        self.infile = ''

    def _set_grids(self, atminifile, ocegridfile):
        """Create grid description files for both atmosphere and ocean"""

        if atminifile:
            self.ATMGRID = self.tempstore.newFile()
            griddes = self.cdo.griddes(input=str(atminifile))
            with open(self.ATMGRID, 'w', encoding='utf-8') as f:
                for line in griddes:
                    print(line, file=f)

        if ocegridfile:
            self.OCEGRID = self.tempstore.newFile()
            griddes = self.cdo.griddes(input=str(ocegridfile))
            with open(self.OCEGRID, 'w', encoding='utf-8') as f:
                for line in griddes:
                    print(line, file=f)

    def _set_atm_fixgrid(self, component, atminifile):
        """Define the command require for correcting model grid"""
        # this could improved using the modelname variable: if EC-Earth, do this...
        if component == 'oifs':
            self.atmfix = f'-setgridtype,regular -setgrid,{self.ATMGRID}'
            self.ATMGRIDAREA = self.cdo.gridarea(input=f'{self.atmfix} {atminifile}')
        elif component == 'cmoratm':
            self.atmfix = ''
            self.ATMGRIDAREA = self.cdo.gridarea(input=f'{atminifile}')
        else:
            sys.exit('Atmospheric component not supported')

    def _set_oce_fixgrid(self, component, ocegridfile, oceareafile):
        """Define the command require for correcting model grid"""

        if ocegridfile and oceareafile:
            # this could improved using the modelname variable: if EC-Earth, do this...
            if component == 'nemo':
                self.OCEGRIDAREA = self.cdo.expr('area=e1t*e2t', input=oceareafile)
                self.ocefix = f'-setgridarea,{self.OCEGRIDAREA}'
            elif component == 'cmoroce':
                self.ocefix = ''
                self.OCEGRIDAREA = self.cdo.gridarea(input=f'{oceareafile}')
            else:
                sys.exit('Oceanic component not supported')

    def set_gridfixes(self, atminifile, ocegridfile, oceareafile, atmcomp, ocecomp):
        """Create all internal grid files and set fixes for atm and oce grids"""
        self._set_grids(atminifile, ocegridfile)
        self._set_atm_fixgrid(atmcomp, atminifile)
        self._set_oce_fixgrid(ocecomp, ocegridfile, oceareafile)

    def make_atm_masks(self, component, atminifile, extra=''):
        """Create land-sea masks for atmosphere model"""
        # prepare ATM LSM: this need to be improved, since it is clearly model dependent
        if component == 'oifs':
            self.LANDMASK = self.cdo.selname('LSM',
                                input=f'-gec,0.5 {extra} {self.atmfix} {atminifile}',
                                options='-t ecmwf -f nc')
            self.SEAMASK = self.cdo.mulc('-1', input=f'-subc,1 {self.LANDMASK}')
        elif component == 'cmoratm':
            self.LANDMASK = self.cdo.selname('sftlf',
                                input=f'-gec,0.5 {extra} {self.atmfix} {atminifile}')
            self.SEAMASK = self.cdo.mulc('-1', input=f'-subc,1 {self.LANDMASK}')

    def make_atm_remap_weights(self, atminifile, remap_method, target):
        """Create atmosphere remap weights"""

        # this creates the method to be called in the remap
        genweight = remap_method.replace('remap', 'gen')
        themethod = getattr(self.cdo(), genweight)

        # self.ATMWEIGHTS = self.cdo.genbil(target, input=f'{self.atmfix} {atminifile}')
        self.ATMWEIGHTS = themethod(target, input=f'{self.atmfix} {atminifile}')
        logging.debug("Atmosphere is remapping with " + genweight + " method")

    def make_oce_remap_weights(self, ocegridfile, remap_method, target):
        """Create ocean remap weights, use remapbil as fallback"""
        if ocegridfile:
            genweight = remap_method.replace('remap', 'gen')
            try:
                themethod = getattr(self.cdo(), genweight)
                self.OCEWEIGHTS = themethod(target, input=f'{self.ocefix} {ocegridfile}')
                logging.debug("Ocean is remapping with " + genweight + " method")
            except:
                themethod = getattr(self.cdo(), 'genbil')
                self.OCEWEIGHTS = themethod(target, input=f'{self.ocefix} {ocegridfile}')
                logging.warning("Ocean is remapping with genbil method cause cannot do " + genweight)

    def chain(self, cmd):
        """Adds a generic cdo operator"""
        self.pipe = '-' + cmd + ' ' + self.pipe

    def setdomain(self, domain):
        """Specify variable domain: used to set needed grid manipulations"""
        self.domain = domain

    def start(self):
        """Cleans pipe for a new application"""
        self.pipe = '{infile}'

    def fixgrid(self, domain=''):
        """Applies grid fixes, requires specifying the domain (atm or oce)"""
        # ocean variables require specifying grid areas
        # atm variables require fixing the grid

        if not domain:
            domain = self.domain

        if not domain:
            sys.exit('You have to define a domain with the setdomain() method first')

        # this should be replaced for a more general "ocean" or "atmosphere"
        if domain == 'oce':
            self.pipe = self.ocefix + ' ' + self.pipe
        elif domain == 'atm':
            self.pipe = self.atmfix + ' ' + self.pipe

    def mask(self, mask_type):
        if mask_type == 'land':
            if not self.LANDMASK:
                sys.exit('Needed grid file not defined, call make_grids method first')
            self.ifthen(self.LANDMASK)
        elif mask_type in ['sea', 'ocean']:
            self.ifthen(self.SEAMASK)

    def masked_meansum(self, mask_type):
        if mask_type == 'land':
            if not self.LANDMASK:
                sys.exit('Needed grid file not defined, call make_grids method first')
            self.chain(f'fldsum -mul {self.ATMGRIDAREA} -mul {self.LANDMASK}')
        elif mask_type in ['sea', 'ocean']:
            self.chain(f'fldsum -mul {self.ATMGRIDAREA} -mul {self.SEAMASK}')
        else:
            self.chain('fldmean')

    def convert(self, offset, factor):
        if offset != 0:
            self.chain(f'addc,{offset}')
        if factor != 1:
            self.chain(f'mulc,{factor}')

    def setname(self, var):
        self.chain(f'setname,{var}')

    def selname(self, var):
        self.chain(f'selname,{var}')

    def selectname(self, var):
        self.pipe = f'-select,name={var} [ ' + self.pipe + ' ] '

    def cat(self):
        self.pipe = '-cat [ ' + self.pipe + ' ] '

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

    def ifthen(self, fname):
        self.pipe = '-ifthen ' + fname + ' ' + self.pipe

    def levels(self, infile):
        return list(map(float, self.cdo.showlevel(input=infile)[0].split()))

    def set_infile(self, infile):
        self.infile = infile

    def execute(self, cmd, *args, input='', keep=False, **kwargs):
        fn = getattr(self.cdo, cmd)
        if not input:
            input = self.infile
        if isinstance(input, list):
            logging.debug('Applying cat: ', input)
            input = '-cat [ ' + ' '.join(input) + ' ]'
        logging.debug('called cdop.execute with: %s', self.pipe)
        logging.debug('input file: %s', input)
        out = fn(input=self.pipe.format(infile=input), *args, **kwargs)
        if not keep:
            self.start()  # clear pipe
        return out

    def output(self, infile, **kwargs):
        logging.debug('call cdop.output with: %s', self.pipe)
        return float(self.execute('output', input=infile, **kwargs)[0])
