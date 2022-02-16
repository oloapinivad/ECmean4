#!/bin/bash

set -ex

# Update the original ECmean post2x2.sh script
# Prepares regridded 2x2 degree files for use with the RK PI scripts
# starting from original ECE4 XIOS data
# The output nc files are in nc4 data

#Updated by P. Davini (ISAC-CNR) - <p.davini@isac.cnr.it> 
#February 2022

if [ $# -ne 3 ]
then
    echo "Usage:   ./post_grid.sh exp YEARSTART YEAREND"
    echo "Example: ./post_grid.sh io01 1990 2000"
    exit 1
fi

# experiment name
expname=$1
# years to be processed
year1=$2
year2=$3

# source config, clean and move to folder
. config_ecmean4.sh
cd $TMPDIR
rm -f $CLIMDIR/*_${expname}_m*_2x2.nc

# create grid descriptor file
$cdo griddes $DATADIR/../../ICMGG${expname}INIT > griddes.txt

# get land sea mask
$cdonc setctomiss,0 -ltc,0.5 -selcode,172 -setgridtype,regular $DATADIR/../../ICMGG${expname}INIT $CLIMDIR/oceanmask_${expname}_2x2.nc

# loop on years
for (( year=$year1; year<=$year2; year++)); do
	
	# extract and interpolate variables
	filein=$DATADIR/${expname}_atm_cmip6_1m_${year}-${year}.nc
	for var in $varlist ; do
		$cdo -L -s cat -setgridtype,regular -setgrid,griddes.txt -selname,$var $filein $CLIMDIR/${var}_${expname}_mon_2x2.nc
	done
done

# mean variables
for var in $varlist ; do
	 [ -f $CLIMDIR/${var}_${expname}_mon_2x2.nc ] && $cdonc timmean $CLIMDIR/${var}_${expname}_mon_2x2.nc $CLIMDIR/${var}_${expname}_mean_2x2.nc
done
	
# back and clean
cd -
rm -rf $TMPDIR
