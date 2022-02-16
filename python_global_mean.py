#!/usr/bin/python3

# this is a tentative python script to convert ECmean global mean operation to python3
# uses a reference file from yaml and cdo bindings (not very efficient)

import yaml
from tabulate import tabulate
from cdo import *
import sys, os, glob
cdo = Cdo()

# arguments
expname = sys.argv[1]
year1 = int(sys.argv[2])
year2 = int(sys.argv[3])

# options (to be moved in a config file?)
ECEDIR=os.path.join("/lus/h2resw01/scratch/ccpd/ece4", expname)
TABDIR=os.path.join(ECEDIR, "ECmean4", "table")
TMPDIR=os.path.join(ECEDIR, "ECmean4", "tmp")
os.makedirs(TABDIR,exist_ok=True)
os.makedirs(TMPDIR,exist_ok=True)
eceinitfile=os.path.join(ECEDIR, "ICMGG" + expname + "INIT")
#cdo.forceOutput = True
#cdo.debug = True

# loop 
varlist=["tas", "psl", "pr", "evspsbl"]
for var in varlist : 
    for year in range(year1,year2+1) :
        print(var, year)
        infile=os.path.join(ECEDIR,"output/oifs", expname + "_atm_cmip6_1m_" + str(year) +  "-" + str(year) +  ".nc")
        outfile=os.path.join(TMPDIR, "tmp_" + str(year) + ".nc")
        cdo.fldmean(input="-setgridtype,regular -setgrid," + eceinitfile + " -selname," +  var + " " + infile, output=outfile)
    cdo.timmean(input="-cat " + os.path.join(TMPDIR, "tmp*.nc"), output=os.path.join(TMPDIR, var + "_mean.nc"))

# reference data
filename="reference.yml"
with open(filename, 'r') as file:
    yummy = yaml.load(file, Loader=yaml.FullLoader)

head=["Var", "Longname", "Units", "ECE4", "OBS", "Obs Dataset"]
cdoformat="%9.4f,1"
global_table=list()

# loop on the variables
for var in varlist :
    filein=os.path.join(TMPDIR, var + "_mean.nc")
    print(filein)
    if os.path.isfile(filein) :
        beta=yummy[var]
        beta['value']=cdo.outputf(cdoformat, input=" -mulc," + beta['factor'] + " " + filein)
        out_sequence=[var , beta['varname'], beta['units'] , beta['value'][0], beta['observations']['val'], beta['observations']['data']]
        global_table.append(out_sequence)

# write the file  with tabulate
tablefile=os.path.join(TABDIR, "Global_Mean_" + expname + "_" + str(year1) + "_" + str(year2) +  ".txt")
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers = head, tablefmt="orgtbl"))

# clean
for f in glob.glob(os.path.join(TMPDIR,"*.nc")):
    os.remove(f)
            
