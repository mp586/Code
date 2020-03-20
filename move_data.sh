#########################################################################
# File Name: move_data.sh
# Author: Wanying Kang
# Created Time: Mon 06 Aug 2018 09:55:08 AM MDT
#########################################################################


for exps in northland_different_albedo northland_same_albedo northland_same_albedo_different_roughness
do
	echo $exps
	cd $exps
	for runs in run0*
	do
		echo $runs
		id="${runs//[!0-9]/}"
		echo $id
	#	cp $runs/atmos_daily.nc ./atmos_daily_$id.nc
		cp $runs/atmos_monthly.nc ./atmos_monthly_$id.nc
	done
	cd -
done
