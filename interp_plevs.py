from isca.util import interpolate_output

# dires_ctl = [
#     'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_commitd15c267',
#     'aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267']

# dires_pert = [
#     'aquaplanet_frierson_insolation_0qflux_mld20_plus_2xCO2_spinup_361_commitd15c267',
#     'square_South_America_frierson_insolation_lepref1_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267',
#     'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_plus_2xCO2_spinup_361_commitd15c267',
#     'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267']
dires_ctl=['square_South_America_lepref07_fixedSSTs_from_realworld_zonallysymm_commit7bb4387', 'square_Africa_lepref07_fixedSSTs_from_realworld_zonallysymm_commit7bb4387', 'two_continents_lepref07_fixedSSTs_from_realworld_zonallysymm_commit7bb4387']
dires_pert = ['square_South_America_lepref07_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387', 'square_Africa_lepref07_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387', 'two_continents_lepref07_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387']

for dire in dires_ctl:
    for i in range(121, 481):
        #print(i)
            infile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly.nc' % i   
            outfile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly_interp.nc' % i
            interpolate_output(infile, outfile, p_levs='INPUT', var_names = ['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum','sphum_u','sphum_v','sphum_w']) #,'ucomp_temp','vcomp_temp','omega_temp','ucomp_height','vcomp_height','omega_height','div'])
    print('Done interpolating '+dire)
for dire in dires_pert:
    for i in range(120, 480):
        #print(i)
            infile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly.nc' % i   
            outfile = '/scratch/mp586/Isca_DATA/ISCA_HPC/'+dire+'/run%04d/atmos_monthly_interp.nc' % i
            interpolate_output(infile, outfile, p_levs='INPUT', var_names = ['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum','sphum_u','sphum_v','sphum_w']) #,'ucomp_temp','vcomp_temp','omega_temp','ucomp_height','vcomp_height','omega_height','div'])
    print('Done interpolating '+dire)
