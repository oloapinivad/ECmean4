#!/bin/bash

DATADIR=/lus/h2resw01/scratch/ccpd/ece4/$expname/output/oifs
CLIMDIR=/lus/h2resw01/scratch/ccpd/ece4/$expname/ecmean4/clim
TMPDIR=/lus/h2resw01/scratch/ccpd/ece4/$expname/ecmean4/tmp_$RANDOM
TABLEDIR=/lus/h2resw01/scratch/ccpd/ece4/$expname/ecmean4/table
mkdir -p $CLIMDIR $TMPDIR $TABLEDIR
cdo=/usr/local/apps/cdo/1.9.10/bin/cdo
cdozip="$cdo -f nc4c -z zip"
cdonc="$cdo -f nc"
remap="remapcon2"
refgrid="r180x90"
varlist="tas psl pr evspsbl clh clm cll rsnt rlnt"
cdoformat="%9.4f,1"
