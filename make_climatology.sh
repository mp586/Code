spin=10
simulyears=30
nyear=$((simulyears-spin))  # total 10 years, minus 2 because the first two years are moved to spinup
# vars=(precipitation t_surf temp ucomp vcomp )
for exps in all_ocean_control
do
	echo $exps
	cd $exps
### move spinup years, need to edit if want to move more than 2 years.
 	mkdir spinup
    spinmon=$((spin*12))
    spinmon=$(printf "%04d" $spinmon)
    for ispin in $(seq -w 0001 $spinmon)
    do
        mv run${ispin} spinup/
    done

### monthly climatology ###
 	ncrcat -O run*/atmos_monthly.nc monthly.nc
 	for mon in `seq 1 12`
 	do
 		echo mon=$mon, monthly climatology
 		monstrt=$((mon-1))
 		mon2=`printf "%02d" $mon`
 		ncra -O -d time,$monstrt,,12 monthly.nc monthlyclimo_${mon2}.nc
 	done
    ncrcat -O monthlyclimo_*.nc monthlyclimo.nc
    cp monthlyclimo.nc /scratch/mp586/Isca_DATA/monthlyclimo_$exps.nc
### end monthly climatology ###

### daily climatology ###
# 	for mon in `seq 1 12`
# 	do
# 		echo mon=$mon, daily climatology
# 		strtrun=`printf "%04d" $((24+mon))`
# 		mon2=`printf "%02d" $mon`
# 		nces -O -n $nyear,4,12 atmos_daily_$strtrun.nc dailyclimo_${mon2}.nc
# 	done
#
#	rm -f dailyclimo.nc dailyclimo_lonavg.nc
#	ncrcat -O dailyclimo_*.nc dailyclimo.nc
#	ncwa -a lon -O dailyclimo.nc dailyclimo_lonavg.nc
#	cp dailyclimo_lonavg.nc /glade/scratch/wykang/dailyclimo_lonavg_$exps.nc
### end daily climatology ###
	cd -
done
