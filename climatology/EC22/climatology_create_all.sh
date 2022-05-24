#!/bin/bash
set -e

# very trivial bash script to create mean and interannual variance for simplified version of reichler and kim built on EC-Mean4
# it includes the variance ratio safecheck
# it requires the data from the required datasets, it is included for reproducibility
# P.Davini, CNR-ISAC, Apr 2022


# to set: time period (default, can be shorter if data are missing)
year1=1990
year2=2019
vars="tas precipitation mean_sea_level_pressure specific_humidity v_component_of_wind u_component_of_wind temperature sosaline sowaflup sozotaux sometauy net_sfc snow_cover analysed_sst sea_ice_fraction"

# cdo details
cdo=cdo205
cdozip="$cdo -f nc4 -z zip"

# targets resolution and years
grids="r180x90 r360x180" # original"
years=${year1}-${year2}


for climtype in EC22 ; do

	# tmpdir
	TMPDIR=$SCRATCH/ecmean_datasets/tmp_$RANDOM
	mkdir -p $TMPDIR

	# dir where ECmean is: define where the climatology will be created
	ECMEANDIR=$SCRATCH/ecmean_datasets/ECmean4

	# loop on vars
	for var in $vars ; do

		echo "Processing $var !!!"

		# var properties (variance ratio to avoid irrelistic values is used)
		case $var in 
			tas) 		dataset=CRU; 	remap_method=remapbil ;;
			precipitation) 	dataset=MSWEP ; remap_method=remapcon ;;
			so*)		dataset=ORAS5 ; remap_method=remapbil ;;
			ileadfra)    	dataset=ORAS5 ; remap_method=remapbil ;;
			specific*)	dataset=ERA5 ;  remap_method=remapbil ;;
			*component*) 	dataset=ERA5 ;  remap_method=remapbil ;;
			temperature) 	dataset=ERA5 ;  remap_method=remapbil ;;
			mean_sea_lev*)	dataset=ERA5 ;  remap_method=remapbil ;;
			snow_cover)	dataset=ERA5 ;  remap_method=remapbil ;;
			sfc_net_tot_all_mon) 	dataset=CERES-EBAF ; remap_method=remapbil ;;
			net_sfc) 	dataset=NOCS ; remap_method=remapbil  ;;
			analysed_sst)	dataset=ESA-CCI-L4 ; remap_method=remapbil  ;;
			sea_ice_fraction)	dataset=ESA-CCI-L4 ; remap_method=remapbil ;;

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

		# safe check for cat on the number of files
		if [[ $(ls $search | wc -l) == 1 ]] ; then
			catcommand=""
		else 
			catcommand="-cat"
		fi

		# cat, select years and yearmean
		tmpfile=$TMPDIR/tmp_${var}.nc 
		rm -f $tmpfile
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
		# check on the number of levels so to implement zonal mean
		cleanvarname=${varname}
		if [[ $levels == 1 ]] ; then
			zoncommand=""
		else 
			zoncommand="-zonmean"
			varname=${varname}_zonal
		fi

	# safecheck on variance: DEPRECATED
	# added condition to set a minimum variance depending on variance ratio, otherwise exclude
	# compute dataset variance and its maximum
	#$cdozip $zoncommand -timvar $varcommand $tmpfile $TMPDIR/variance_tmp_${var}.nc
	#maxvar=$(cdo output -fldmax -vertmax $TMPDIR/variance_tmp_${var}.nc)
	# define minimum accepted variance as a variance_ratio of the max variance 
	#factor=$(awk -v a="$maxvar" -v b="${variance_ratio}" 'BEGIN{print (a / b)}')
	#echo "Minimum accepted variance: $factor"

		# create grid files
		for grid in $grids ; do 
			mkdir -p $grid
		
	        	if [[ $grid == "original" ]] ; then
        	        	remap=""
       			else
        	        	remap="-${remap_method},$grid"
        		fi

			# safecheck on variance, version 2
			# use 5 sigma from the mean of the log10 distribution
			# special treatment of weird values for SST
	        	if [[ $var == "analysed_sst" ]] ; then
	                	cleaning="-setrtomiss,0,5e-03"
	        	else
	                	cleaning="-setctomiss,0"
	        	fi

			# compute the variance on the needed grid, and move it to log10
			$cdozip -b 64 log10 ${cleaning} -timvar ${zoncommand} ${remap} ${varcommand} ${tmpfile} $TMPDIR/variance_tmp_${var}.nc

			# mean field, without area or level wieghts 
	        	meanvar=$(cdo output -fldmean,weights=FALSE -vertmean,weights=FALSE  $TMPDIR/variance_tmp_${var}.nc)

			# more complicated for std, use a false area and expr to estimate it also for vertical profiles
			cdo addc,1 -mulc,0 $TMPDIR/variance_tmp_${var}.nc $TMPDIR/false_area.nc
			sdvar=$(cdo output -expr,"out=sqrt(fldmean(vertmean(sqr($cleanvarname-fldmean(vertmean($cleanvarname))))))" -setgridarea,$TMPDIR/false_area.nc $TMPDIR/variance_tmp_${var}.nc)
       			# define minimum accepted variance as a variance_ratio of the max variance
        		echo $meanvar $sdvar

			# extract range
        		lowfactor=$(awk -v a="$meanvar" -v b="${sdvar}" 'BEGIN{print (10^(a - 5*b))}')
        		highfactor=$(awk -v a="$meanvar" -v b="${sdvar}" 'BEGIN{print (10^(a + 5*b))}')
        		echo $lowfactor $highfactor


			# loop on climate and variance files
			for kind in climate variance ; do

				if [[ $kind == climate ]] ; then
					timcommand="-timmean"
					factorcommand=""
				elif [[ $kind == variance ]] ; then
					timcommand="-timvar"
					if [[ $climtype == "EC22_nofilter" ]] ; then
						factorcommand=""
					elif [[ $climtype == "EC22" ]] ; then
						factorcommand="-setvrange,$lowfactor,$highfactor"
						#factorcommand="-setrtomiss,0,$factor"
					fi
				fi
				newsuffix=${varname}_${dataset}_${grid}_${firstyear}-${lastyear}.nc
				CLIMDIR=$ECMEANDIR/climatology/${climtype}/${grid}
                        	mkdir -p $CLIMDIR
				echo $cdozip ${factorcommand} ${cleaning} ${timcommand} ${zoncommand} ${remap} ${varcommand} $tmpfile  $CLIMDIR/${kind}_${newsuffix}
				$cdozip ${factorcommand} ${cleaning} ${timcommand} ${zoncommand} ${remap} ${varcommand} $tmpfile  $CLIMDIR/${kind}_${newsuffix}

			done
		done
	# clean
	rm $tmpfile 
	done
rm -rf $TMPDIR
done
