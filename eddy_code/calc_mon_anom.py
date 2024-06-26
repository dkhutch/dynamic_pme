#!/usr/bin/env python
import numpy as np
import netCDF4 as nc
import os

infile = 'plevel.nc'
outfile = 'anom.nc'
monfile = 'mon_av.nc'

indirs = []
fnames = []
for i in range(380,390):
    indirs.append('output%3d' % i)
    fnames.append('output%3d/month_plevel.nc' % i)

ndirs = len(indirs)
year_per_dir = 5
nyrs = ndirs * year_per_dir

fnames_join = ' '.join(fnames)

if not os.path.exists(monfile):
    print('Making monthly averages')
    cmd = f'ncrcat -O {fnames_join} mon_cat.nc'
    os.system(cmd) 
    cmd = f'av_mon.py mon_cat.nc {monfile}'
    os.system(cmd)
    os.remove('mon_cat.nc')
    print('Made monthly averages')
else:
    print('Monthly average already exists')

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

print('loaded monthly data')

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
    fo.history = 'calc_mon_anom.py'
    
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

    u_o = fo.createVariable('u_anom', 'f4', ('time', 'level', 'lat', 'lon'))
    u_o.units = 'm/s'
    u_o.long_name = 'zonal wind'

    v_o = fo.createVariable('v_anom', 'f4', ('time', 'level', 'lat', 'lon'))
    v_o.units = 'm/s'
    v_o.long_name = 'meridional wind'

    q_o = fo.createVariable('q_anom', 'f4', ('time', 'level', 'lat', 'lon'))
    q_o.units = 'kg/kg'
    q_o.long_name = 'specific humidity'

    for yr in range(nyr):
        print('year', yr+1)
        yrind = yr * 365
        for t in range(nmon):
            print('month', t+1)
            t0 = tind[t] + yrind
            t1 = tind[t+1] + yrind
            q = f.variables['sphum'][t0:t1,:,:,:]
            u = f.variables['ucomp'][t0:t1,:,:,:]
            v = f.variables['vcomp'][t0:t1,:,:,:]
            tt = f.variables['time'][t0:t1]

            qanom = q - q_mon[t,:,:,:]
            uanom = u - u_mon[t,:,:,:]
            vanom = v - v_mon[t,:,:,:]

            q_o[t0:t1,:,:,:] = qanom[:]
            u_o[t0:t1,:,:,:] = uanom[:]
            v_o[t0:t1,:,:,:] = vanom[:]
            time_o[t0:t1] = tt[:]

    fo.close()
    f.close()
    os.chdir('..')