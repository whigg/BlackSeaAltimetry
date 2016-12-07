#!/usr/bin/env python
#
# blacksea_plot_meanEKE.py
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
resdir='/home/ctroupin/DIVA/BlackSea_Altimetry/results/V1/Filtered/EKE/'
resdir='/home/ctroupin/DIVA/BlackSea_Altimetry/results/V2_seminorm/Unfiltered/EKE/'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/EKE/V2_seminorm/Unfiltered/'
coastdir='/home/ctroupin/DataOceano/Coastlines/'
coastfile='blacksea_coast2.dat'
resbasename='mean_EKE_'
figtype='.png'

yearmin=2010
yearmax=2013

# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=27.,42,40.,47.,0.,.0
vmin,vmax,dvar,dvar2=0.,100.,2.5,25.

resolution=300

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
cmap=plt.cm.hot_r

# Norm and levels to plot
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.001)
levels2plot=np.arange(vmin-0.000,vmax+0.000,dvar)
bounds2= np.arange(vmin-0.000,vmax+0.000,dvar2)
filelist=sorted(glob.glob(resdir+str(yearmin)+'/'+resbasename+'*'))

lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
lonc[lonc<-10]=np.nan
latc[latc<-10]=np.nan
xc,yc=m(lonc,latc)

# --------------------------------
# Load coordinates 
# --------------------------------

firstfile=resdir+resbasename+str(yearmin)+'.nc'
nc=netcdf.Dataset(firstfile)
lon1 = nc.variables['longitude'][:]
lat1 = nc.variables['latitude'][:]
nc.close()


# Projection
llon1, llat1  = np.meshgrid(lon1,lat1)
lon1,lat1= m(llon1,llat1)

#
EKEmean = np.zeros_like(llon1)
yearlist=range(yearmin,yearmax)
nyears=len(yearlist)
for years in yearlist:

    file2load=resdir+resbasename+str(years)+'.nc'
    print os.path.basename(file2load)
    
    with netcdf.Dataset(file2load) as nc:
        EKE=nc.variables['EKE'][:].squeeze()
        
    EKEmean+=EKE/nyears
        
    # --------------
    # Make the plots
    #---------------

    fig=plt.figure()
            

    ax = fig.add_subplot()
    plt.title('EKE'+' ('+str(years)+')')
    contour1=m.contourf(lon1,lat1,EKE,levels2plot,cmap=cmap,norm=norm,extend='max')
    m.fillcontinents(color='w',zorder=3)
    m.plot(xc,yc,'k',lw=0.5,zorder=4)
    cbar=fig.colorbar(contour1,orientation='vertical',shrink=0.7,pad=0.04,aspect=15)
    cbar.set_label('($cm^2/s^2$)',fontname='Times New Roman',fontsize=18,rotation=0)
    cbar.set_ticks(bounds2)
    cbar.ax.set_xticklabels(bounds2,fontname='Times New Roman',fontsize=16)

    m.drawparallels(np.arange(latmin-1,latmax+1,2), linewidth=0.5,zorder=2,
                            labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16)
    m.drawmeridians(np.arange(lonmin,lonmax+0.0001,2), linewidth=0.5,zorder=2,
                            labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16)

    plt.savefig(figdir+'EKE'+str(years)+figtype, dpi=resolution, facecolor='w', edgecolor='w',
       transparent=True, bbox_inches='tight', pad_inches=0.1)
    ##        
    #plt.show()
    plt.close()
##
##
fig=plt.figure()
            

ax = fig.add_subplot()
plt.title('EKE'+' ('+str(yearmin)+'-'+str(yearmax)+')')
contour1=m.contourf(lon1,lat1,EKEmean,levels2plot,cmap=cmap,norm=norm,extend='max')
m.fillcontinents(color='w',zorder=3)
m.plot(xc,yc,'k',lw=0.5,zorder=4)
cbar=fig.colorbar(contour1,orientation='vertical',shrink=0.7,pad=0.04,aspect=15)
cbar.set_label('($cm^2/s^2$)',fontname='Times New Roman',fontsize=18,rotation=0)
cbar.set_ticks(bounds2)
cbar.ax.set_xticklabels(bounds2,fontname='Times New Roman',fontsize=16)

m.drawparallels(np.arange(latmin-1,latmax+1,2), linewidth=0.5,zorder=2,
                        labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16)
m.drawmeridians(np.arange(lonmin,lonmax+0.0001,2), linewidth=0.5,zorder=2,
                        labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16)

plt.savefig(figdir+'EKE'+str(yearmin)+'_'+str(yearmax)+figtype, dpi=resolution, facecolor='w', edgecolor='w',
   transparent=True, bbox_inches='tight', pad_inches=0.1)
##        
#plt.show()
plt.close()
