from isca.util import interpolate_output

dires_ctl = ['narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_six_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387']

dires_pert = ['narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]

# dires_ctl=['square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387']
# dires_pert=['square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387']

for dire in dires_ctl:
    for i in range(121, 481): #range(121, 481) # range(481, 541)
        #print(i)
            # infile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_daily.nc' % i   
            # outfile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_daily_interp.nc' % i
            # interpolate_output(infile, outfile, p_levs='INPUT', var_names = ['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum','sphum_u','sphum_v','sphum_w']) #,'ucomp_temp','vcomp_temp','omega_temp','ucomp_height','vcomp_height','omega_height','div'])
            infile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly.nc' % i   
            outfile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly_interp.nc' % i
            interpolate_output(infile, outfile, p_levs='INPUT', var_names = ['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum','sphum_u','sphum_v','sphum_w','div','vor']) #,'ucomp_temp','vcomp_temp','omega_temp','ucomp_height','vcomp_height','omega_height','div'])
    print('Done interpolating '+dire)
for dire in dires_pert:
    for i in range(120,480):#range(120, 480) # range(481, 541)
        #print(i)
            # infile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_daily.nc' % i   
            # outfile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_daily_interp.nc' % i
            # interpolate_output(infile, outfile, p_levs='INPUT', var_names = ['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum','sphum_u','sphum_v','sphum_w']) #,'ucomp_temp','vcomp_temp','omega_temp','ucomp_height','vcomp_height','omega_height','div'])
            infile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly.nc' % i   
            outfile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly_interp.nc' % i
            interpolate_output(infile, outfile, p_levs='INPUT', var_names = ['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum','sphum_u','sphum_v','sphum_w','div','vor']) #,'ucomp_temp','vcomp_temp','omega_temp','ucomp_height','vcomp_height','omega_height','div'])
    print('Done interpolating '+dire)
