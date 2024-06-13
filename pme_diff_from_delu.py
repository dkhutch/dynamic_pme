#!/usr/bin/env python
import numpy as np
import netCDF4 as nc
import argparse

parser = argparse.ArgumentParser(description='calculate P-E using flux divergence')
parser.add_argument('control', help='control experiment file containing ucomp, vcomp, sphum on pressure levels')
parser.add_argument('pert', help='perturbation experiment file containing ucomp, vcomp, sphum on pressure levels')
parser.add_argument('output', help='output file with P-E derived field')
args = parser.parse_args()

edgefile = 'latlon_edge.nc'

f = nc.Dataset(args.control,'r')
q = f.variables['sphum'][:]
u_ctl = f.variables['ucomp'][:]
v_ctl = f.variables['vcomp'][:]
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]
pres = f.variables['level'][:]
f.close()

f = nc.Dataset(args.pert, 'r')
u_pert = f.variables['ucomp'][:]
v_pert = f.variables['vcomp'][:]
f.close()

pres = pres * 100 # convert to Pa

nlat, nlon = lat.shape[0], lon.shape[0]

f = nc.Dataset(edgefile,'r')
latb = f.variables['latb'][:]
lonb = f.variables['lonb'][:]
f.close()

rad_earth = 6376000
grav = 9.81
# rho_w = 997

uq_int = np.trapz((u_pert - u_ctl) * q, x=pres, axis=1)
vq_int = np.trapz((v_pert - v_ctl) * q, x=pres, axis=1) 

lonmat, latmat = np.meshgrid(lon, lat)
coslat = np.cos(np.deg2rad(latmat))

duq_x = np.gradient(uq_int, axis=2)
dvq_y = np.gradient(vq_int, axis=1)

dlat = np.deg2rad(np.diff(latb))
dlon = np.deg2rad(np.diff(lonb))

dx = dlon * coslat * rad_earth
dy = dlat * rad_earth

dx = dx[np.newaxis, :, :]
dy = dy[np.newaxis, :, np.newaxis]

grad_x = duq_x / dx
grad_y = dvq_y / dy

pme_calc = (grad_x + grad_y) / grav

f = nc.Dataset(args.output, 'w')
f.history = f'pme_diff_from_delu.py {args.control} {args.pert} {args.output}'

f.createDimension('time', 0)
f.createDimension('lat', nlat)
f.createDimension('lon', nlon)

lat_o = f.createVariable('lat', 'f8', ('lat'))
lat_o.units = 'degrees_north'
lat_o[:] = lat[:]

lon_o = f.createVariable('lon', 'f8', ('lon'))
lon_o.units = 'degrees_east'
lon_o[:] = lon[:]

pme_o = f.createVariable('pme_calc', 'f8', ('time','lat','lon'))
pme_o.units = 'kg/m2/s'
pme_o[:] = pme_calc[:]

f.close()