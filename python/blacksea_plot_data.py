#!/usr/bin/env python

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors

#-------------------
# Plot the DIVA data
#-------------------

os.system('clear')
basemap_resolution = 'c'


#filedir='/home/ctroupin/DataOceano/AVISO/BlackSea/4-time_weight_10.0/'
filedir='/home/ctroupin/DataOceano/AVISO/BlackSea/2-mission_concat/'
#filebasename='dt_blacksea_merged_sla_vxxc_201210*'
filebasename='dt_blacksea_*sla_vxxc_201210*.dat'

coastdir='/home/ctroupin/DataOceano/Coastlines/'
coastfile='blacksea_coast2.dat'
#figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/data/gathered/'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/data/daily/'
figtype='.png'

# Select colormaps
cmap=plt.cm.spectral

# Region of interest
lonmin,lonmax,latmin,latmax=27,42.,40.,48.

# Spacing between x/y labels
dlon,dlat = 3.0,2.0

# Limits of the colorbar
vmin,vmax,dvar = -20.0,20.0,5.0
bounds2 = np.arange(vmin,vmax+0.0001,dvar)
norm = colors.Normalize(vmin=vmin,vmax=vmax)


#------------------------------------------------------------------------------

# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

# Create basemap then used for the plot
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            lat_ts=0.5*(latmin+latmax),\
            resolution=basemap_resolution)

# Load extracted coastline
valexc=-999.
lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
lonc=np.ma.masked_where(lonc==valexc,lonc)
latc=np.ma.masked_where(latc==valexc,latc)
xc,yc=m(lonc,latc)

# Start the loop on the data files
# present in the directory

for infile in sorted(glob.glob(filedir + filebasename)):

    filename= os.path.basename(infile)
    print('Working on file ' + filename)
    figname=figdir+filename[0:-4]+figtype

    # Load data from file
    file2load=filedir+filename


    # Check if file is not empty
    if(os.path.getsize(file2load) > 0):
        # Here we only read the first three columns to
        # avoid reading strings corresponding to satellites
        #
        
        lon,lat,field,=np.loadtxt(file2load,usecols=(0,1,2),unpack='True')
        field*=100.        # convert meters in centimeters

        # Make the plot
        fig=plt.figure()
        ax = fig.add_axes([0.1,0.2,0.8,0.7])


        m.ax=ax
        x,y = m(lon, lat)
        scat = m.scatter(x, y, c=field,s=5,zorder=2, edgecolor='none',norm=norm,cmap=cmap)
        m.plot(xc,yc,lw=0.5,color='k')
##        m.drawcoastlines(ax=ax)
##        m.fillcontinents(color='black', ax=ax)

        m.drawparallels(np.arange(latmin,latmax,dlat), linewidth=0,labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(lonmin,lonmax,dlon), linewidth=0,labels=[0, 0, 1, 0])

        #m.bluemarble()
        # Add the colorbar

        cbar=fig.colorbar(scat,cmap=cmap,orientation='horizontal',shrink=0.75,pad=0.02,extend='both')
        cbar.set_label('(cm)',fontsize=18)
        cbar.set_ticks(bounds2)
        cbar.ax.set_xticklabels(bounds2,fontsize=16)

        plt.savefig(figname, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        plt.close()
    else:
        print 'file ',filename,' is empty'
        
