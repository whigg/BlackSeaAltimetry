#!/usr/bin/env python
#
# blacksea_check_time.py
#
#

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


filedir='/home/ctroupin/DataOceano/AVISO/BlackSea/4-time_weight_10.0/'
filebasename='dt_blacksea_merged_sla_vxxc_1993070*'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/data/'
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
        
        lon,lat,field,weight,time=np.loadtxt(file2load,usecols=(0,1,2,3,5),unpack='True')
        field*=100.        # convert meters in centimeters

        # Make the plot
        fig=plt.figure()
        ax = fig.add_axes([0.1,0.2,0.8,0.7])
    
        plt.plot(time-time.min(),weight)
        plt.show()
        plt.close()
    else:
        print 'file ',filename,' is empty'
        
