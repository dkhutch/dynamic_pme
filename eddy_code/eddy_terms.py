#!/usr/bin/env python
#PBS -P y99
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l storage=gdata/hh5+scratch/y99+scratch/v45
#PBS -l wd
#PBS -j oe
#PBS -N eddy_prism

import numpy as np
import netCDF4 as nc
import os

infile = 'anom.nc'
outfile = 'eddy_term.nc'

indirs = []
fnames = []
for i in range(380,390):
    indirs.append('output%3d' % i)
    fnames.append('output%3d/anom.nc' % i)

ndirs = len(indirs)
year_per_dir = 5
nyrs = ndirs * year_per_dir

monfile = 'mon_av.nc'

f = nc.Dataset(monfile,'r')
q_mon = f.variables['sphum'][:]
u_mon = f.variables['ucomp'][:]
v_mon = f.variables['vcomp'][:]
pres = f.variables['level'][:]
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
latb = f.variables['latb'][:]
lonb = f.variables['lonb'][:]
f.close()

nlev, nlat, nlon = pres.shape[0], lat.shape[0], lon.shape[0]
nmon = 12
months = np.array([0,31,28,31,30,31,30,31,31,30,31,30,31], 'i4') # include zero here for indexing
tind = np.cumsum(months)
nyr = 5

for i in range(ndirs):
    os.chdir(indirs[i])
    print('dir', indirs[i])
    f = nc.Dataset(infile,'r')
    fo = nc.Dataset(outfile,'w')
    fo.history = 'eddy_terms.py'

    fo.createDimension('time',0)
    fo.createDimension('level', nlev)
    fo.createDimension('lat', nlat)
    fo.createDimension('lon', nlon)
    fo.createDimension('latb', nlat+1)
    fo.createDimension('lonb', nlon+1)

    time_o = fo.createVariable('time', 'f8', ('time'))
    time_o.units = 'days since 0001-01-01 00:00:00'
    time_o.calendar = 'NOLEAP'

    pres_o = fo.createVariable('level', 'f8', ('level'))
    pres_o.units = 'hPa'
    pres_o[:] = pres[:]

    lat_o = fo.createVariable('lat', 'f8', ('lat'))
    lat_o.units = 'degrees_north'
    lat_o[:] = lat[:]

    lon_o = fo.createVariable('lon', 'f8', ('lon'))
    lon_o.units = 'degrees_east'
    lon_o[:] = lon[:]

    latb_o = fo.createVariable('latb', 'f8', ('latb'))
    latb_o.units = 'degrees_north'
    latb_o[:] = latb[:]

    lonb_o = fo.createVariable('lonb', 'f8', ('lonb'))
    lonb_o.units = 'degrees_east'
    lonb_o[:] = lonb[:]

    u_o = fo.createVariable('uqbar', 'f4', ('time', 'level', 'lat', 'lon'))
    u_o.units = 'm/s'
    u_o.long_name = 'zonal wind anomaly * specific humidity anomaly'

    v_o = fo.createVariable('vqbar', 'f4', ('time', 'level', 'lat', 'lon'))
    v_o.units = 'm/s'
    v_o.long_name = 'meridional wind anomaly * specific humidity anomaly'

    for yr in range(nyr):
        print('year', yr+1)
        yrind = yr * 365
        for t in range(nmon):
            t_out_ind = yr * nmon + t
            print('month', t+1)
            t0 = tind[t] + yrind
            t1 = tind[t+1] + yrind
            q = f.variables['q_anom'][t0:t1,:,:,:]
            u = f.variables['u_anom'][t0:t1,:,:,:]
            v = f.variables['v_anom'][t0:t1,:,:,:]
            tt = f.variables['time'][t0:t1]

            uqbar = np.ma.mean(u * q, axis=0)
            vqbar = np.ma.mean(v * q, axis=0)
            tbar = np.ma.mean(tt, axis=0)

            time_o[t_out_ind] = tt
            u_o[t_out_ind, :, :, :] = uqbar
            v_o[t_out_ind, :, :, :] = vqbar
    
    f.close()
    fo.close()
    os.chdir('..')


