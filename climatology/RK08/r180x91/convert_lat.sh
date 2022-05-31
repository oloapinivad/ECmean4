#!/bin/bash

filelist="*.nc"
for file in $filelist ; do 
	check=$(cdo griddes $file | grep yinc | cut -f2 -d"=")
	echo $check
	if [[ $((check)) == -2 ]] ; then
		echo "need to convert $file"
		cdo -f nc4 -z zip invertlat $file tmp.nc 
		mv tmp.nc $file
	else 
		echo "File is ok $file"
		#cdo -f nc4 -z zip copy $file tmp.nc
		#mv tmp.nc $file
	fi
	
done

