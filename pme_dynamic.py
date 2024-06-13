import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util as cutl
import cmaps

plio_file = '../P_minus_E/pme_eplio.nc'
glac_file = '../P_minus_E/pme_glac.nc'
dyn_pme_file = 'pme_diff_delu_glac.nc'

outfile_jja = 'pme_dynamic_jja.pdf'
outfile_djf = 'pme_dynamic_djf.pdf'

f = nc.Dataset(plio_file, 'r')
evap_plio = f.variables['evap'][:] 
precip_plio = f.variables['precip'][:] 
temp_plio = f.variables['t_ref'][:] 
lat = f.variables['lat'][:] 
lon = f.variables['lon'][:] 
f.close()

f = nc.Dataset(glac_file, 'r') 
evap_glac = f.variables['evap'][:]
precip_glac = f.variables['precip'][:] 
temp_glac = f.variables['t_ref'][:] 
f.close()

f = nc.Dataset(dyn_pme_file,'r')
pme_dyn = f.variables['pme_calc'][:]
f.close()

evap_plio = evap_plio 
precip_plio = precip_plio 
evap_glac = evap_glac 
precip_glac = precip_glac 

pme_plio = precip_plio - evap_plio
pme_glac = precip_glac - evap_glac

delta_pme = pme_glac - pme_plio 
delta_T = temp_glac - temp_plio

delta_T_av = np.mean(delta_T)

alpha = 0.07

delta_held = alpha * delta_T * pme_plio

jja = [5,6,7]
djf = [0,1,11]

delta_held_jja = np.mean(delta_held[jja,:,:], axis=0) * 86400.
delta_held_djf = np.mean(delta_held[djf,:,:], axis=0) * 86400.
delta_pme_jja = np.mean(delta_pme[jja,:,:], axis=0) * 86400.
delta_pme_djf = np.mean(delta_pme[djf,:,:], axis=0) * 86400.
pme_plio_jja = np.mean(pme_plio[jja,:,:], axis=0) * 86400.
pme_plio_djf = np.mean(pme_plio[djf,:,:], axis=0) * 86400.
pme_dyn_jja = np.mean(pme_dyn[jja,:,:], axis=0) * 86400.
pme_dyn_djf = np.mean(pme_dyn[djf,:,:], axis=0) * 86400.

delta_held_jja_z = np.mean(delta_held_jja, axis=1)
delta_held_djf_z = np.mean(delta_held_djf, axis=1)
delta_pme_jja_z = np.mean(delta_pme_jja, axis=1)
delta_pme_djf_z = np.mean(delta_pme_djf, axis=1)
pme_plio_jja_z = np.mean(pme_plio_jja, axis=1)
pme_plio_djf_z = np.mean(pme_plio_djf, axis=1)

delta_held_jja, lon1 = cutl.add_cyclic_point(delta_held_jja, coord=lon) 
delta_held_djf, _ = cutl.add_cyclic_point(delta_held_djf, coord=lon) 
delta_pme_jja, _ = cutl.add_cyclic_point(delta_pme_jja, coord=lon) 
delta_pme_djf, _ = cutl.add_cyclic_point(delta_pme_djf, coord=lon)
pme_plio_jja, _ = cutl.add_cyclic_point(pme_plio_jja, coord=lon)
pme_plio_djf, _ = cutl.add_cyclic_point(pme_plio_djf, coord=lon) 
pme_dyn_jja, _ = cutl.add_cyclic_point(pme_dyn_jja, coord=lon)
pme_dyn_djf, _ = cutl.add_cyclic_point(pme_dyn_djf, coord=lon)

amin = -1.
amax = 1.
ivl = .1
levels = np.arange(amin, amax+ivl/2., ivl)

bmin = -10.
bmax = 10.
bivl = 1.
blevels = np.arange(bmin, bmax+bivl/2., bivl)

fig = plt.figure(num=1, figsize=(10,6))


ax1 = fig.add_subplot(2,2,1, projection=ccrs.PlateCarree())

# can also use cmap=cmaps.cmocean_balance
# h1 = ax1.contourf(lon1, lat, pme_plio_jja, blevels, extend="both", cmap=cmaps.BlueWhiteOrangeRed)
h1 = ax1.contourf(lon1, lat, pme_plio_jja, blevels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h1, fraction=0.046, pad=0.04)
ax1.set_title('Early Plio P-E JJA (mm/d)')
ax1.coastlines()
gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False


ax1 = fig.add_subplot(2,2,2, projection=ccrs.PlateCarree())

# can also use cmap=cmaps.cmocean_balance
h1 = ax1.contourf(lon1, lat, delta_held_jja, levels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h1, fraction=0.046, pad=0.04)
ax1.set_title('Held-Soden 06 $\Delta$(P-E) JJA (mm/d)')
ax1.coastlines()
gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

ax1 = fig.add_subplot(2,2,3, projection=ccrs.PlateCarree())

h1 = ax1.contourf(lon1, lat, pme_dyn_jja, levels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h1, fraction=0.046, pad=0.04)
ax1.set_title('Dynamic $\Delta$(P-E) JJA (mm/d)')
ax1.coastlines()
gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False


ax2 = fig.add_subplot(2,2,4, projection=ccrs.PlateCarree())
h2 = ax2.contourf(lon1, lat, delta_pme_jja, levels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h2, fraction=0.046, pad=0.04)
ax2.set_title('Total $\Delta$(P-E) JJA (mm/d)')
ax2.coastlines()
gl = ax2.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

plt.savefig(outfile_jja)










fig = plt.figure(num=2, figsize=(10,6))


ax1 = fig.add_subplot(2,2,1, projection=ccrs.PlateCarree())

h1 = ax1.contourf(lon1, lat, pme_plio_djf, blevels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h1, fraction=0.046, pad=0.04)
ax1.set_title('Early Plio P-E DJF (mm/d)')
ax1.coastlines()
gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False


ax1 = fig.add_subplot(2,2,2, projection=ccrs.PlateCarree())

h1 = ax1.contourf(lon1, lat, delta_held_djf, levels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h1, fraction=0.046, pad=0.04)
ax1.set_title('Held-Soden 06 $\Delta$(P-E) DJF (mm/d)')
ax1.coastlines()
gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

ax1 = fig.add_subplot(2,2,3, projection=ccrs.PlateCarree())

h1 = ax1.contourf(lon1, lat, pme_dyn_djf, levels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h1, fraction=0.046, pad=0.04)
ax1.set_title('Dynamic $\Delta$(P-E) DJF (mm/d)')
ax1.coastlines()
gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

ax2 = fig.add_subplot(2,2,4, projection=ccrs.PlateCarree())
h2 = ax2.contourf(lon1, lat, delta_pme_djf, levels, extend="both", cmap=cmaps.MPL_BrBG)
cbar = plt.colorbar(h2, fraction=0.046, pad=0.04)
ax2.set_title('Total $\Delta$(P-E) DJF (mm/d)')
ax2.coastlines()
gl = ax2.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

plt.savefig(outfile_djf)

