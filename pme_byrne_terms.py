#!/usr/bin/env python
import numpy as np
import netCDF4 as nc

'''
Calculate the 4 terms needed for equation A5 of 
Byrne and O'Gorman 2015.
http://dx.doi.org/10.1175/JCLI-D-15-0369.s1
'''

ctl_av_file = 'uvq_av_eplio.nc'
pert_av_file = 'uvq_av_glac.nc'
ctl_eddy_file = 'eddy_av_eplio.nc'
pert_eddy_file = 'eddy_av_glac.nc'
outfile = 'pme_byrne_glac_minus_eplio.nc'


f = nc.Dataset(ctl_av_file,'r')
q_ctl = f.variables['sphum'][:]
u_ctl = f.variables['ucomp'][:]
v_ctl = f.variables['vcomp'][:]
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]
pres = f.variables['level'][:]
latb = f.variables['latb'][:]
lonb = f.variables['lonb'][:]
f.close()

f = nc.Dataset(pert_av_file,'r')
q_pert = f.variables['sphum'][:]
u_pert = f.variables['ucomp'][:]
v_pert = f.variables['vcomp'][:]
f.close()

f = nc.Dataset(ctl_eddy_file, 'r')
uqbar_ctl = f.variables['uqbar'][:]
vqbar_ctl = f.variables['vqbar'][:]
f.close()

f = nc.Dataset(pert_eddy_file, 'r')
uqbar_pert = f.variables['uqbar'][:]
vqbar_pert = f.variables['vqbar'][:]
f.close()

#---------------------------------------
# Calculate grid terms for differentials
#---------------------------------------
pres = pres * 100 # convert to Pa
nlat, nlon = lat.shape[0], lon.shape[0]

rad_earth = 6376000
grav = 9.8 # GFDL CM2.1 uses exactly 9.8 for gravity
lonmat, latmat = np.meshgrid(lon, lat)
coslat = np.cos(np.deg2rad(latmat))

dlat = np.deg2rad(np.diff(latb))
dlon = np.deg2rad(np.diff(lonb))

dx = dlon * coslat * rad_earth
dy = dlat * rad_earth

dx = dx[np.newaxis, :, :]
dy = dy[np.newaxis, :, np.newaxis]

#---------------------------------------------
# 1. Calculate mean thermodynamic change δ MTh 
#---------------------------------------------

uq_int = np.trapz(u_ctl * (q_pert - q_ctl), x=pres, axis=1)
vq_int = np.trapz(v_ctl * (q_pert - q_ctl), x=pres, axis=1) 
duq_x = np.gradient(uq_int, axis=2)
dvq_y = np.gradient(vq_int, axis=1)
grad_x = duq_x / dx
grad_y = dvq_y / dy

pme_MTh = (grad_x + grad_y) / grav

#-----------------------------------------
# 2. Calculate mean dynamic change δ MDyn 
#-----------------------------------------

uq_int = np.trapz((u_pert - u_ctl) * q_ctl, x=pres, axis=1)
vq_int = np.trapz((v_pert - v_ctl) * q_ctl, x=pres, axis=1) 
duq_x = np.gradient(uq_int, axis=2)
dvq_y = np.gradient(vq_int, axis=1)
grad_x = duq_x / dx
grad_y = dvq_y / dy

pme_MDyn = (grad_x + grad_y) / grav

#-------------------------------
# 3. Calculate eddy term δ Eddy 
#-------------------------------

uq_ctl = np.trapz(uqbar_ctl, x=pres, axis=1)
vq_ctl = np.trapz(vqbar_ctl, x=pres, axis=1) 
uq_pert = np.trapz(uqbar_pert, x=pres, axis=1)
vq_pert = np.trapz(vqbar_pert, x=pres, axis=1)
duq = uq_pert - uq_ctl
dvq = vq_pert - vq_ctl
duq_x = np.gradient(duq, axis=2)
dvq_y = np.gradient(dvq, axis=1)
grad_x = duq_x / dx
grad_y = dvq_y / dy

pme_Eddy = (grad_x + grad_y) / grav

#-------------------------------
# 4. Calculate Non-Linear term
#-------------------------------

uq_int = np.trapz((u_pert - u_ctl) * (q_pert - q_ctl), x=pres, axis=1)
vq_int = np.trapz((v_pert - v_ctl) * (q_pert - q_ctl), x=pres, axis=1) 
duq_x = np.gradient(uq_int, axis=2)
dvq_y = np.gradient(vq_int, axis=1)
grad_x = duq_x / dx
grad_y = dvq_y / dy

pme_NL = (grad_x + grad_y) / grav

#-------------------------------
# Write output
#-------------------------------

f = nc.Dataset(outfile, 'w')
f.history = f'pme_byrne_terms.py'

f.createDimension('time', 0)
f.createDimension('lat', nlat)
f.createDimension('lon', nlon)

time_o = f.createVariable('time', 'f8', ('time'))
time_o.units = 'months'
time_o[:] = np.arange(1,13)

lat_o = f.createVariable('lat', 'f8', ('lat'))
lat_o.units = 'degrees_north'
lat_o[:] = lat[:]

lon_o = f.createVariable('lon', 'f8', ('lon'))
lon_o.units = 'degrees_east'
lon_o[:] = lon[:]

pme_1 = f.createVariable('pme_mth', 'f8', ('time','lat','lon'), fill_value=-1.e10)
pme_1.units = 'kg/m2/s'
pme_1.long_name = 'Mean Thermodynamic term of d(P-E)'
pme_1[:] = pme_MTh[:]

pme_2 = f.createVariable('pme_mdyn', 'f8', ('time','lat','lon'), fill_value=-1.e10)
pme_2.units = 'kg/m2/s'
pme_2.long_name = 'Mean Dynamic term of d(P-E)'
pme_2[:] = pme_MDyn[:]

pme_3 = f.createVariable('pme_eddy', 'f8', ('time','lat','lon'), fill_value=-1.e10)
pme_3.units = 'kg/m2/s'
pme_3.long_name = 'Eddy term of d(P-E)'
pme_3[:] = pme_Eddy[:]

pme_4 = f.createVariable('pme_nl', 'f8', ('time','lat','lon'), fill_value=-1.e10)
pme_4.units = 'kg/m2/s'
pme_4.long_name = 'Non-Linear term of d(P-E)'
pme_4[:] = pme_NL[:]


f.close()