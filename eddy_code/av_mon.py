#!/usr/bin/env python

import numpy as np 
import netCDF4 as nc 
import sys
import argparse

parser = argparse.ArgumentParser(description="Convert monthly series to monthly climatology")
parser.add_argument("origin", help="input: multi-year monthly average netCDF file")
parser.add_argument("destination", help="output: climatological monthly average netCDF file")

args = parser.parse_args()

ncfile_o = nc.Dataset(args.origin, 'r')
ncfile_d = nc.Dataset(args.destination, 'w')

glob_atts = ncfile_o.__dict__
newhist = 'av_mon.py {} {} \n '.format(args.origin, args.destination) 
if 'history' in glob_atts:
    glob_atts['history'] = newhist + glob_atts['history']
else:
    glob_atts['history'] = newhist
ncfile_d.setncatts(glob_atts)

unlimdimname = False
unlimdim = None

for dimname,dim in ncfile_o.dimensions.items():
    if dim.isunlimited():
        unlimdimname = dimname
        unlimdim = dim
        nt = len(unlimdim)
        nyr = int(nt / 12)
        ncfile_d.createDimension(dimname,12)
    else:
        ncfile_d.createDimension(dimname,len(dim))

months = np.array([31,28,31,30,31,30,31,31,30,31,30,31], 'f8')

assert unlimdim is not None, 'Time dimension must be a record dimension (unlimited)'

varnames = ncfile_o.variables.keys()

for varname in varnames:
    ncvar = ncfile_o.variables[varname]
    data = ncvar[:] 

    if hasattr(ncvar, '_FillValue'):
        FillValue = ncvar._FillValue
    else:
        FillValue = None

    if unlimdimname in ncvar.dimensions:
        if len(data.shape) > 1:
            keepdims = data.shape[1:]
            data = np.reshape(data, list(np.hstack(([nyr, 12], keepdims[:]))) )
        else: # then it's a 1D time variable
            data = np.reshape(data, [nyr, 12])
        data = np.ma.mean(data, axis=0)
    
    var = ncfile_d.createVariable(varname, ncvar.dtype, ncvar.dimensions, fill_value=FillValue,
          zlib=True, complevel=5)

    attdict = ncvar.__dict__ 
    # remove the fillvalue attribute and others we don't want
    for rematt in ['_FillValue','time_avg_info','missing_value','cell_methods','standard_name']:
        if rematt in attdict: 
            del attdict[rematt]
    var.setncatts(attdict)

    var[:] = data[:]    
