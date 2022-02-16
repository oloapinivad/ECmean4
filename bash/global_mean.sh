#!/bin/bash

set -ex

# Derived from original ECmean script
# Computes global averages of radiative fluxes (plus a few selected fields) 
# February 2022

if [ $# -ne 1 ]
then
    echo "Usage:   ./global_mean.sh EXP"
    echo "Example: ./global_mean.sh io01"
    exit 1
fi

# experiment name
expname=$1
. config_ecmean4.sh


#year1=$(cdo showyear -seltimestep,1 $CLIMDIR/tas_${expname}_mean_2x2.nc)
tablename=${TABLEDIR}/Global_Mean_Table_${expname}.txt
rm -f $tablename


net_toa=$( $cdo outputf,$cdoformat -fldmean -add $CLIMDIR/rsnt_${expname}_mean_2x2.nc $CLIMDIR/rlnt_${expname}_mean_2x2.nc )

varrad="rlnt rsnt net_toa"
# radiation values
echo -e "Radiation" >> ${tablename}
echo -e "Variable\tunits\t\tECE4-${expname}\t\tTrenberth09" >> $tablename

for var in $varrad ; do

	filein=$CLIMDIR/${var}_${expname}_mean_2x2.nc
	
	case $var in 
	
		"rsnt")		varname="TOA net SW" ;          trenval=239.4 ;         expval=$( $cdo outputf,$cdoformat -fldmean $filein ) ;;
        	"rlnt")	 	varname="TOA net LW" ;          trenval=-238.5 ;        expval=$( $cdo outputf,$cdoformat -fldmean $filein ) ;;
		"net_toa")      varname="Net TOA   " ; 		trenval=0.9 ;		expval=${net_toa} ;;

	esac

	echo -e "${varname}\tW/m2\t\t${expval}\t\t${trenval}" >> $tablename

done

# global mean values
echo -e >> $tablename 
echo -e "Global Mean" >> ${tablename}
echo -e "Variable\tunits\t\tECE4-${expname}\t\tObservations" >> $tablename

varmean="tas psl pr evspsbl cll clm clh"
for var in $varmean ; do

	eraval=""
    	datas=""
	case $var in

        	"tas")          varname="Air T at 2m    ";      units="K         ";	oper="" ;    	eraval=287.575 ; 	datas="ERAI(1990-2010)";;
		"psl")          varname="MSLP           ";      units="Pa        ";     oper="" ;	eraval=101135  ; 	datas="ERAI(1990-2010)";;
		"pr")		varname="Tot Precipit.  ";      units="mm/day    ";    	oper="-mulc,86400" ; 	eraval=2.92339 ; 	datas="ERAI(1990-2010)";;
		"evspsbl")      varname="Evaporation    ";      units="mm/day    ";	oper="-mulc,86400" ;;    
		"cll")          varname="Low CC         ";      units="0-1       ";    	oper="-divc,100" ;	eraval=0.36627 ; 	datas="ERAI(1979-2006)";;
        	"clm")          varname="Medium CC      ";      units="0-1       ";     oper="-divc,100" ; 	eraval=0.17451 ; 	datas="ERAI(1979-2006)";;
        	"clh")          varname="High CC        ";      units="0-1       ";     oper="-divc,100" ; 	eraval=0.28723 ; 	datas="ERAI(1979-2006)";;
		

	esac
	expval=$($cdo outputf,$cdoformat $oper -fldmean $CLIMDIR/${var}_${expname}_mean_2x2.nc)
	echo -e "${varname}\t${units}\t${expval}\t\t${eraval}\t$datas" >> $tablename
done


