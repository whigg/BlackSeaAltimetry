#!/usr/bin/env python
#
# blacksea_plot_data_results.py
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
resdir='/home/ctroupin/DIVA/BlackSea_Altimetry/results/test_param/'
datadir='/home/ctroupin/DataOceano/AVISO/BlackSea/Filtered/3-time_concat_60days/'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/results/test_param/'
#resbasename='dt_blacksea_'
resbasename='results'
figtype='.png'

yearmin=2012
yearmax=2012

# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=27.,42,40.,47.,0.,.0
vmin,vmax,dvar=-3.,3.,0.5

resolution=150

# Create basemap then used for the plot
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            lat_ts=0.5*(latmin+latmax),\
            resolution=basemap_resolution)

# ----------------------------------------------------------------------------

# Create directory if necessary
if not os.path.exists(figdir):
    os.makedirs(figdir)
        
# Colormap
cmap=plt.cm.RdYlBu_r

# Norm and levels to plot
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.001)
levels2plot=np.arange(vmin-0.000,vmax+0.000,dvar)


# --------------------------------
# Load coordinates 
# --------------------------------

#firstfile=sorted(glob.glob(resdir+str(yearmin)+'/'+resbasename+'*'))[0]
firstfile=sorted(glob.glob(resdir+resbasename+'*'))[0]
nc=netcdf.Dataset(firstfile)
lon1 = nc.variables['x'][:]
lat1 = nc.variables['y'][:]
nc.close()


# Projection
llon1, llat1  = np.meshgrid(lon1,lat1)
lon1,lat1= m(llon1,llat1)

for years in range(yearmin,yearmax+1):
    filelist1 = sorted(glob.glob(resdir+str(years)+'/'+resbasename+'*'))
    
    for i in range(0,10):#,len(filelist1)):

        print i

        # Load file content
        file2load1=filelist1[i]
        figtime=os.path.basename(file2load1).split('_')[5]
        nc=netcdf.Dataset(file2load1)
        sla1=nc.variables['analyzed_field'][:]*100.
        nc.close()  

        #datafile=glob.glob(datadir+'dt_blacksea_'+'*'+figtime+'*')[0]
        lond,latd,data=np.loadtxt(datafile,usecols=(0,1,2),unpack=True)
        #lond,latd=m(lond,latd)                                                                

        
        # Compute gradient
        [sla1_x,sla1_y]=np.gradient(sla1)
        
        # --------------
        # Make the plots
        #---------------

        fig=plt.figure(num=None, figsize=(10,5), facecolor='w', edgecolor='k')
                
        
        ax = fig.add_subplot(121)
        plt.title('x-gradient')
        contour1=m.contourf(lon1,lat1,sla1_x,levels2plot,cmap=cmap,norm=norm,extend='both')
        #m.plot(lond,latd,'ko',ms=0.5)
        plt.colorbar(contour1,orientation='horizontal')
        
        ax = fig.add_subplot(122)
        plt.title('y-gradient')
        contour2=m.contourf(lon1,lat1,sla1_y,levels2plot,cmap=cmap,norm=norm,extend='both')
        #m.plot(lond,latd,'ko',ms=0.5)
        plt.colorbar(contour2,orientation='horizontal')
        plt.savefig(figdir+resbasename+figtime+figtype, dpi=resolution, facecolor='w', edgecolor='w',
          transparent=True, bbox_inches='tight', pad_inches=0.)
        
        #plt.show()
        plt.close()


