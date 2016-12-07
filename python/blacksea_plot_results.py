#!/usr/bin/env python
#
# blacksea_plot_results.py
#
#
# ---------------------------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from mpl_toolkits.basemap import Basemap
from matplotlib import rc
from matplotlib  import colors
import time, calendar
from datetime import date
import sys

t0 = time.time()

basemap_resolution = 'l'
resdir='/home/ctroupin/DIVA/BlackSea_Altimetry/results/V2/'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/results/V2/'
avisodir='/home/ctroupin/DataOceano/AVISO/BlackSea/NetCDF_gridded/'
resbasename='dt_blacksea_merged_sla_vxxc_2011'
avisobasename='dt_blacksea_allsat_msla_h_'
figtype='.png'

yearmin=2011
yearmax=2011

# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=27.,42,40.,47.,0.,.0
vmin,vmax,dvar=-10.,10.,0.5

resolution=150


# ----------------------------------------------------------------------------

# Create directory if necessary
if not os.path.exists(figdir):
    os.makedirs(figdir)
        
# Initialise basemap projection
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            lat_ts=0.5*(latmin+latmax),\
            resolution=basemap_resolution)

# Colormap
cmap=plt.cm.RdYlBu_r

# Norm and levels to plot
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.001)
levels2plot=np.arange(vmin-0.000,vmax+0.000,dvar)

# --------------------------------
# Load coordinates 
# --------------------------------

firstfile=sorted(glob.glob(resdir+str(yearmin)+'/'+resbasename+'*'))[0]
nc=netcdf.Dataset(firstfile)
lon1 = nc.variables['x'][:]
lat1 = nc.variables['y'][:]
nc.close()

firstfile=sorted(glob.glob(avisodir+str(yearmin)+'/'+avisobasename+'*'))[0]
nc=netcdf.Dataset(firstfile)
lon2 = nc.variables['lon'][:]
lat2 = nc.variables['lat'][:]
nc.close()

# Projection
llon1, llat1  = np.meshgrid(lon1,lat1)
lon1,lat1 = m(llon1, llat1)
llon2, llat2  = np.meshgrid(lon2,lat2)
lon2,lat2 = m(llon2, llat2)

for years in range(yearmin,yearmax+1):
    filelist1 = sorted(glob.glob(resdir+str(years)+'/'+resbasename+'*'))
    filelist2 = sorted(glob.glob(avisodir+str(years)+'/'+avisobasename+'*'))
    
    for i in range(0,1):#len(filelist1)-30):

        print i

        # Load file content
        file2load1=filelist1[i]
        nc=netcdf.Dataset(file2load1)
        sla1=nc.variables['analyzed_field'][:]*100.
        nc.close()

        file2load2=filelist2[i+30]
        nc=netcdf.Dataset(file2load2)
        sla2=np.squeeze(nc.variables['sla'][:])*100.
        nc.close()
       
        # --------------
        # Make the plots
        #---------------

        fig=plt.figure()
        fig.patch.set_facecolor('white')
        fig.patch.set_alpha(0.)
        plt.title('sdvndfsnvfn')
        ax = fig.add_subplot(121)
        m.ax=ax
        contour1=m.contourf(lon1,lat1,sla1-sla1.mean(),levels2plot,cmap=cmap,norm=norm,extend='both')

        ax = fig.add_subplot(122)
        m.ax=ax
        contour2=m.contourf(lon2,lat2,sla2-sla2.mean(),levels2plot,cmap=cmap,norm=norm,extend='both')

##        plt.savefig(figdir+figname++'_'+figtime+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##          transparent=True, bbox_inches='tight', pad_inches=0.)
        plt.colorbar(contour1,orientation='horizontal')
        plt.show()
        plt.close()


