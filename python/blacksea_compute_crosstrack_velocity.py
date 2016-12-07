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
from blacksea_functions import calc_uv_track

t0 = time.time()

basemap_resolution = 'l'
coastdir='/home/ctroupin/DataOceano/Coastlines/'
coastfile='blacksea_coast2.dat'
datadir='/home/ctroupin/DataOceano/AVISO/BlackSea/Filtered/2-mission_concat/'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/crosstrack_velocity/'
databasename='dt_blacksea_'
figtype='.png'
figbasename='cross_track_velocity_vfec_'

vmin,vmax=0.0,0.5
cmap=plt.cm.hot_r
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.001)

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

# coastline
valexc=-999.
lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
lonc=np.ma.masked_where(lonc==valexc,lonc)
latc=np.ma.masked_where(latc==valexc,latc)

# ----------------------------------------------------------------------------

# Create directory if necessary
if not os.path.exists(figdir):
    os.makedirs(figdir)
        
for years in range(yearmin,yearmax+1):

    for month in range(1,13):
        mmonth=str(month).zfill(2)
        filelist1 = sorted(glob.glob(datadir+databasename+'*'+'_'+str(years)+mmonth+'*'))

        # Create empty arrays for lon and lat
        lonvel,latvel,vel=np.array([]),np.array([])
        
        fig=plt.figure()
        ax = fig.add_subplot(111)
            
        for i in range(0,len(filelist1)):
            file2load=filelist1[i]

            #print 'Working on '+os.path.basename(file2load)
            num_lines = sum(1 for line in open(file2load))

            
            
            if num_lines>10:
                lon,lat,sla,track,cycle=np.loadtxt(file2load,usecols=(0,1,2,6,7),unpack=True)
                track_u=np.unique(track)
                ntracks=len(track_u)

                for t in range(0,ntracks):
                    #print 'Working on track '+ str(track_u[t])
                    goodtrack=np.where(track==track_u[t])
                    geostr_vel,ugeostr,vgeostr,distSLA,alpha=calc_uv_track(lon[goodtrack],lat[goodtrack],sla[goodtrack])
                    
                    # scatter plot
                    lon2plot=0.5*(lon[goodtrack][:-1]+lon[goodtrack][1:])
                    lat2plot=0.5*(lat[goodtrack][:-1]+lat[goodtrack][1:])
                    scat1=plt.scatter(lon2plot,lat2plot,s=5,c=abs(geostr_vel),edgecolor='none',cmap=cmap,norm=norm)

                    # create vectors to later store data
                    lonvel=np.concatenate((lonvel,lon))
                    latvel=np.concatenate((latvel,lat))
                    vel=np.concatenate((vel,geostr_vel))
                    
        cbar=plt.colorbar()
        cbar.set_label('m/s',rotation=0)
        plt.plot(lonc,latc,'k',lw=0.5)
        ax.set_xlim(27.,42.)
        ax.set_ylim(40.,48.)
        plt.title(calendar.month_name[month]+' '+str(years),fontsize=20)
        
                  
        figname=figdir+figbasename+str(years)+mmonth+figtype
        plt.savefig(figname, facecolor='w', edgecolor='w',bbox_inches='tight',pad_inches=0.)
        #plt.show()
        plt.close()

      
