import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util as cutl
import cmaps

infile = 'pme_byrne_glac_minus_eplio.nc'
plot_base = 'pme_4terms_'

f = nc.Dataset(infile,'r')
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
pme_mth = f.variables['pme_mth'][:]
pme_mdyn = f.variables['pme_mdyn'][:]
pme_eddy = f.variables['pme_eddy'][:]
pme_nl = f.variables['pme_nl'][:]
f.close()

months = np.array([31,28,31,30,31,30,31,31,30,31,30,31], 'f8')
jja = [5,6,7]
djf = [0,1,11]
ann = list(range(12))

time_indices = [jja, djf, ann]
time_labels = ['JJA','DJF','ANN']

levels = np.round(np.linspace(-1., 1., 21), decimals=1)

def make_plot_av(var_in, tind):
    var_av = np.ma.average(var_in[tind,:,:], axis=0, weights=months[tind])
    var_av = var_av * 86400.
    var_av, lon1 = cutl.add_cyclic_point(var_av, coord=lon)
    return var_av, lon1

for t in range(3):
    tind = time_indices[t]
    print('calculating means for', time_labels[t])
    pme_mth_av, lon1 = make_plot_av(pme_mth, tind)
    pme_mdyn_av, _ = make_plot_av(pme_mdyn, tind)
    pme_eddy_av, _ = make_plot_av(pme_eddy, tind)
    pme_nl_av, _ = make_plot_av(pme_nl, tind)

    print('making plots for', time_labels[t])

    fig = plt.figure(num=t, figsize=(10,6))

    ax1 = fig.add_subplot(2,2,1, projection=ccrs.PlateCarree())

    h1 = ax1.contourf(lon1, lat, pme_mth_av, levels, extend='both', cmap=cmaps.MPL_BrBG)
    ax1.set_title('Thermodynamic $\Delta$ (P-E) Glac - Early Plio ' + time_labels[t])
    ax1.coastlines()
    gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    ax1 = fig.add_subplot(2,2,2, projection=ccrs.PlateCarree())

    h1 = ax1.contourf(lon1, lat, pme_mdyn_av, levels, extend='both', cmap=cmaps.MPL_BrBG)
    ax1.set_title('Dynamic $\Delta$ (P-E) Glac - Early Plio ' + time_labels[t])
    ax1.coastlines()
    gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    ax1 = fig.add_subplot(2,2,3, projection=ccrs.PlateCarree())

    h1 = ax1.contourf(lon1, lat, pme_eddy_av, levels, extend='both', cmap=cmaps.MPL_BrBG)
    ax1.set_title('Eddy $\Delta$ (P-E) Glac - Early Plio ' + time_labels[t])
    ax1.coastlines()
    gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    ax1 = fig.add_subplot(2,2,4, projection=ccrs.PlateCarree())

    h1 = ax1.contourf(lon1, lat, pme_nl_av, levels, extend='both', cmap=cmaps.MPL_BrBG)
    ax1.set_title('Non-Linear $\Delta$ (P-E) Glac - Early Plio ' + time_labels[t])
    ax1.coastlines()
    gl = ax1.gridlines(draw_labels=True, alpha=0.5, color='black', linewidth=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False    

    fig.subplots_adjust(bottom=0.2)
    cbar1 = fig.add_axes([0.12,0.1,0.78,0.04])
    clb = fig.colorbar(h1, cax=cbar1, orientation='horizontal')
    clb.set_label('mm/day')
    labs = np.array(levels[::2])
    clb.set_ticks(labs)
    clb.set_ticklabels(labs)

    plt.savefig(plot_base + time_labels[t] + '.pdf')

 