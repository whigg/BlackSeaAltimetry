#--------------------------------------------------------------
#
# blacksea_EKE_netcdf.py
#
#
#
# ctroupin, Summer 2014
#--------------------------------------------------------------

import os
import glob 
import numpy as np
import numexpr as ne
import netCDF4 as netcdf
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib  import colors
from scipy import interpolate
from scipy import ndimage
import datetime, calendar
import time as ttime

from blacksea_functions import *

inputdir='/home/ctroupin/DIVA/BlackSea_Altimetry/results/V2_seminorm/Unfiltered/'
yearmin,yearmax=2010,2013
outputdir=inputdir+'/EKE/'
resbasename='dt_blacksea_merged_sla_vxxc_'
sec2day=86400.

if not(os.path.exists(outputdir)):
    os.makedirs(outputdir)

cmap = plt.cm.hot_r
vmin,vmax,dvar=0.,250.,10.
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.001)
levels2plot=np.arange(vmin-0.000,vmax+0.000,dvar)

#---------------------------------------------
# Get grid information from first file in list
# Load SLA
#---------------------------------------------

SLAdir=inputdir+str(yearmin)+'/'
filelist = sorted(glob.glob(SLAdir + resbasename +'*'))
nc=netcdf.Dataset(filelist[0])
lonSLA = nc.variables['x'][:]
latSLA = nc.variables['y'][:]
SLA = nc.variables['analyzed_field'][:]
nc.close()

valex=SLA.fill_value
print 'valex = '+str(valex)

#plt.pcolormesh(SLA);plt.show()

llonSLA,llatSLA=np.meshgrid(lonSLA,latSLA)

#---------------------------
# Create a mask and erode it
#---------------------------

mask=np.ones(SLA.shape)
mask[SLA<=valex]=0
mask=ndimage.binary_erosion(mask)

umask, vmask=uvpmask(mask)

#SLA = np.ma.masked_where(SLA==valex, SLA)


#------------------------------------------
# Get Coriolis
#------------------------------------------
f = np.sin(np.deg2rad(llatSLA))
f *= 4.
f *= np.pi
f /= sec2day
              
lonu = half_interp(llonSLA[:,:-1], llonSLA[:,1:])
latu = half_interp(llatSLA[:,:-1], llatSLA[:,1:])
lonv = half_interp(llonSLA[:-1], llonSLA[1:])
latv = half_interp(llatSLA[:-1], llatSLA[1:])

#------------------------------------------
# Get pm and pn
#------------------------------------------
pm = np.zeros_like(llonSLA)
pm[:,1:-1] = distLonLat(lonu[:,:-1], latu[:,:-1],
                                     lonu[:,1:],  latu[:,1:])
pm[:,0] = pm[:,1]
pm[:,-1] = pm[:,-2]
 
pn = np.zeros_like(llonSLA)
pn[1:-1] = distLonLat(lonv[:-1], latv[:-1],
                                   lonv[1:],  latv[1:])
pn[0] = pn[1]
pn[-1] = pn[-2]
       
np.reciprocal(pm, out=pm)
np.reciprocal(pn, out=pn)



#--------------------------------
# Start the loop on the SLA files
#--------------------------------

# Loop on years
for years in range(yearmin,yearmax):
    SLAdir=inputdir+str(years)+'/'
    filelist = sorted(glob.glob(SLAdir + resbasename +'*'))
    nfiles=len(filelist)

    # Allocate
    EKEmean=np.zeros_like(SLA)
    umean=np.zeros_like(SLA)
    vmean=np.zeros_like(SLA)

    for SLAfile in filelist:
        (a,filename)= os.path.split(SLAfile)
   
        filedate = filename.split('_')[5]
        print filedate
        year=int(filedate[:4])
        month=int(filedate[4:6])
        day=int(filedate[6:])

        filedatesec = calendar.timegm((year,month,day,0,0,0,))/sec2day

        # Load SLA
        nc=netcdf.Dataset(SLAfile)
        SLA = nc.variables['analyzed_field'][:]
        nc.close()
        SLA = np.ma.masked_where(SLA==valex, SLA)
    
        ugv, vgv= getSurfGeostrVel(f, SLA, pm, pn, umask, vmask)
        ugv[np.logical_or(mask==0,abs(ugv)>10)]=valex
        vgv[np.logical_or(mask==0,abs(vgv)>10)]=valex

        # EKE computed using the velocity components derived from SLA
        EKE=ne.evaluate('5000*(ugv**2+vgv**2)')

        EKEmean+=EKE/nfiles
        umean+=ugv/nfiles
        vmean+=vgv/nfiles           
            
    #--------------------
    # Create a new NetCDF
    #--------------------

    EKEmean=np.ma.masked_where(EKEmean>5000.0,EKEmean)
    EKEmean=np.ma.masked_where(np.logical_and(llonSLA>34.0,llatSLA>45.0),EKEmean)

##    fig=plt.figure()
##    ax = fig.add_subplot()
##    plt.title('EKE'+' ('+str(years)+')')
##    contour1=plt.contourf(lonSLA,latSLA,EKEmean,levels2plot,cmap=cmap,norm=norm,extend='max')
##    plt.show(); plt.show()
    
    create_netcdf_eke(outputdir+'mean_EKE_'+str(years)+'.nc',len(latSLA),len(lonSLA),valex,filedatesec,\
                  lonSLA,latSLA,umean,vmean,EKEmean)


print 'Files written in '+outputdir
print ' '



