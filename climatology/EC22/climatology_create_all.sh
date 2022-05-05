#!/bin/bash
set -e

# very trivial bash script to create mean and interannual variance for simplified version of reichler and kim built on EC-Mean4
# it includes the variance ratio safecheck
# it requires the data from the required datasets, it is included for reproducibility
# P.Davini, CNR-ISAC, Apr 2022


# to set: time period (default, can be shorter if data are missing)
year1=1990
year2=2019
vars="mean_sea_level_pressure specific_humidity v_component_of_wind u_component_of_wind temperature sosaline sosstsst ileadfra sowaflup sozotaux sometauy precipitation net_sfc tas"
vars="snow_cover"
vars="analysed_sst sea_ice_fraction"
vars="tas"

# directory where the data to create the climatology are
TMPDIR=/work/scratch/users/paolo/ecmean_datasets/tmp_$RANDOM
ECMEANDIR=/work/scratch/users/paolo/ecmean_datasets/ECmean4
mkdir -p $TMPDIR


# cdo details
cdozip="cdo205 -f nc4 -z zip"

# targets resolution and years
grids="r180x90 r360x180" # original"
years=${year1}-${year2}

# loop on vars
for var in $vars ; do

	echo "Processing $var !!!"

	# var properties (variance ratio to avoid irrelistic values is used)
	case $var in 
		tas) 		dataset=CRU; 	remap_method=remapbil ; variance_ratio=1e3 ;;
		precipitation) 	dataset=MSWEP ; remap_method=remapcon ; variance_ratio=1e6 ;;
		so*)		dataset=ORAS5 ; remap_method=remapbil ; variance_ratio=1e6 ;;
		ileadfra)    	dataset=ORAS5 ; remap_method=remapbil ; variance_ratio=1e6 ;;
		specific*)	dataset=ERA5 ;  remap_method=remapbil ; variance_ratio=1e9 ;;
		*component*) 	dataset=ERA5 ;  remap_method=remapbil ; variance_ratio=1e6 ;;
		temperature) 	dataset=ERA5 ;  remap_method=remapbil ; variance_ratio=1e6 ;;
		mean_sea_lev*)	dataset=ERA5 ;  remap_method=remapbil ; variance_ratio=1e3 ;;
		snow_cover)	dataset=ERA5 ;  remap_method=remapbil ; variance_ratio=1e3 ;;
		sfc_net_tot_all_mon) 	dataset=CERES-EBAF ; remap_method=remapbil ; variance_ratio=1e6 ;;
		net_sfc) 	dataset=NOCS ; remap_method=remapbil ; variance_ratio=1e6 ;;
		analysed_sst)	dataset=ESA-CCI-L4 ; remap_method=remapbil ; variance_ratio=1e3 ;;
		sea_ice_fraction)	dataset=ESA-CCI-L4 ; remap_method=remapbil ; variance_ratio=1e6 ;;

	esac

	# dataset directories: can use the DATADIR but can specify other details
	DATADIR=/work/scratch/users/paolo/ecmean_datasets
	case $dataset in 
		CRU) 	search="$DATADIR/$dataset/data/*.nc" ;;
		MSWEP)	search="$DATADIR/$dataset/data/*.nc" ;;
		ORAS5)	search="$DATADIR/$dataset/data/$var*.nc" ;; 
		ERA5)   search="/work/datasets/obs/ERA5/$var/mon/${dataset}_${var}*.nc" ;; 
		CERES-EBAF) search="$DATADIR/$dataset/data/*.nc" ;;
		NOCS) search="$DATADIR/$dataset/data/*net_sfc*.nc" ;;
		ESA-CCI-L4) search="/work/datasets/obs/$dataset/day/${var}/*${var}*.nc" ;;
	esac

	# clean temp
	rm -f tmp_${var}.nc

	# safe check for cat on the number of files
	if [[ $(ls $search | wc -l) == 1 ]] ; then
		catcommand=""
	else 
		catcommand="-cat"
	fi

	# cat, select years and yearmean
	tmpfile=$TMPDIR/tmp_${var}.nc 
	cdo yearmean -selyear,$year1/$year2 $catcommand $search $tmpfile

	# extract real first and last year
	firstyear=$(cdo showyear $tmpfile | head -n1 | cut -f2 -d" ")
	lastyear=$(cdo showyear $tmpfile | tail -n2 | rev | cut -f1 -d" " | rev)
	echo $firstyear $lastyear


	# specific properties for ERA5
	if [[ $dataset == "ERA5"  ]] ; then
		varname=$(cdo -t ecmwf showname $tmpfile | xargs)
		varcommand="-setname,$varname"
	elif [[ $dataset == "CRU" ]] ; then
		varname=$var
		varcommand="-setname,$varname -selname,tmp"
	else
		varname=$var
		varcommand=""
	fi

	levels=$(cdo nlevel $tmpfile | head -n1)
	# check on the number of levels
	if [[ $levels == 1 ]] ; then
		zoncommand=""
	else 
		zoncommand="zonmean"
		varname=${varname}_zonal
	fi

	# safecheck on variance
	# added condition to set a minimum variance depending on variance ratio, otherwise exclude
	# compute dataset variance and its maximum
	$cdozip $zoncommand -timvar $varcommand $tmpfile $TMPDIR/variance_tmp_${var}.nc
	maxvar=$(cdo output -fldmax -vertmax $TMPDIR/variance_tmp_${var}.nc)
	# define minimum accepted variance as a variance_ratio of the max variance 
	factor=$(awk -v a="$maxvar" -v b="${variance_ratio}" 'BEGIN{print (a / b)}')
	echo "Minimum accepted variance: $factor"

	# create grid files
	for grid in $grids ; do 
		mkdir -p $grid
		
	        if [[ $grid == "original" ]] ; then
        	        remap=""
       		else
        	        remap="-${remap_method},$grid"
        	fi

		# loop on climate and variance files
		for kind in climate variance ; do

			if [[ $kind == climate ]] ; then
				timcommand="-timmean"
				factorcommand=""
			elif [[ $kind == variance ]] ; then
				timcommand="-timvar"
				factorcommand="-setrtomiss,0,$factor"
			fi
			newsuffix=${varname}_${dataset}_${grid}_${firstyear}-${lastyear}.nc
			CLIMDIR=$ECMEANDIR/climatology/EC22/${grid}
                        mkdir -p $CLIMDIR
			$cdozip $zoncommand $factorcommand ${varcommand} $timcommand $remap $tmpfile  $CLIMDIR/${kind}_${newsuffix}

		done
	done
	# clean
	rm $tmpfile 

done
rm -rf $TMPDIR

