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


function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}
