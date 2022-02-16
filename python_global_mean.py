#!/usr/bin/python3

# this is a tentative python script to convert ECmean global mean operation to python3
# uses a reference file from yaml and cdo bindings (not very efficient)

import sys
import os
import glob
import yaml
from tabulate import tabulate
from cdo import *
cdo = Cdo()

# arguments
expname = sys.argv[1]
year1 = int(sys.argv[2])
year2 = int(sys.argv[3])

# options (to be moved in a config file?)
ECEDIR = os.path.join("/lus/h2resw01/scratch/ccpd/ece4", expname)
TABDIR = os.path.join(ECEDIR, "ECmean4", "table")
TMPDIR = os.path.join(ECEDIR, "ECmean4", "tmp")
os.makedirs(TABDIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)
eceinitfile = os.path.join(ECEDIR, "ICMGG" + expname + "INIT")
#cdo.forceOutput = True
#cdo.debug = True

# loop
var_field = ["tas", "psl", "pr", "evspsbl", "cll", "clm", "clh"]
var_radiation = ["rlnt", "rsnt", "rsns", "rlns", "hfss", "hfls"]
for var in var_field + var_radiation:
    for year in range(year1, year2+1):
        print(var, year)
        infile = os.path.join(ECEDIR, "output/oifs", expname +
                              "_atm_cmip6_1m_" + str(year) + "-" + str(year) + ".nc")
        outfile = os.path.join(TMPDIR, "tmp_" + str(year) + ".nc")
        cdo.fldmean(input="-setgridtype,regular -setgrid," +
                    eceinitfile + " -selname," + var + " " + infile, output=outfile)
    cdo.timmean(input="-cat " + os.path.join(TMPDIR, "tmp*.nc"),
                output=os.path.join(TMPDIR, var + "_mean.nc"))

# var_extra: NET_TOA
extra_radiation=["net_toa", "net_sfc"]

# this must be done with python, it is too dumb with CDO
cdo.add(input=os.path.join(TMPDIR, "rsnt" + "_mean.nc ") + os.path.join(TMPDIR,
        "rlnt" + "_mean.nc"), output=os.path.join(TMPDIR, "net_toa" + "_mean.nc"))
cdo.add(input=os.path.join(TMPDIR, "rsns" + "_mean.nc ") + " -add " + os.path.join(TMPDIR,
        "rlns" + "_mean.nc ") + " -add " + os.path.join(TMPDIR,
        "hfss" + "_mean.nc ") + os.path.join(TMPDIR, "hfls" + "_mean.nc "), 
        output=os.path.join(TMPDIR, "net_sfc" + "_mean.nc"))

# reference data: it is badly written but it can be implemented in a much more intelligent
# and modulable way
filename = "reference.yml"
with open(filename, 'r') as file:
    yummy = yaml.load(file, Loader=yaml.FullLoader)

# define options for the output table
head = ["Var", "Longname", "Units", "ECE4", "OBS", "Obs Dataset"]
cdoformat = "%9.4f,1"
global_table = list()
radiation_table = list()

# loop on the variables
for var in var_field + var_radiation + extra_radiation:
    filein = os.path.join(TMPDIR, var + "_mean.nc")
    print(filein)
    if os.path.isfile(filein):
        beta = yummy[var]
        # also this is better with python (do we need iris? I don't thnk so)
        beta['value'] = cdo.outputf(
            cdoformat, input=" -mulc," + beta['factor'] + " " + filein)
        out_sequence = [var, beta['varname'], beta['units'], beta['value']
                        [0], beta['observations']['val'], beta['observations']['data']]
        global_table.append(out_sequence)

# write the file  with tabulate: cool python feature
tablefile = os.path.join(TABDIR, "Global_Mean_" +
                         expname + "_" + str(year1) + "_" + str(year2) + ".txt")
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers=head, tablefmt="orgtbl"))

# clean
for f in glob.glob(os.path.join(TMPDIR, "*.nc")):
    os.remove(f)
