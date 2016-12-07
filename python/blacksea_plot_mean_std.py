#!/usr/bin/env python
#
# compute and plot stats from the results
# obtained with Diva
#
# --------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import calendar

resolution=300
resultbasedir='/data_local/DIVA/BlackSea_Altimetry/results/V1/Filtered/'
resultbasefile='dt_blacksea_merged_sla_vfec_*'
avisobasedir='/data_local/AVISO/BlackSea/NetCDF_gridded/'
avisobasefile='dt_blacksea_allsat_msla_h_*'
yearstart=1993
yearend=2012
figdir='/home/ctroupin/Presentations/MyOcean_ScienceDays/PosterDiva/figures2/'
figname1='blacksea_mean_max_min'
figname2='blacksea_mean_max_min_zoom'

nyears=1+yearend-yearstart
ndays=int(np.ceil(nyears*365.25))

fieldmin_vec=np.zeros(ndays)
fieldmax_vec=np.zeros(ndays)
fieldmean_vec=np.zeros(ndays)
fieldrms_vec=np.zeros(ndays)

fieldmin_vec2=np.zeros(ndays)
fieldmax_vec2=np.zeros(ndays)
fieldmean_vec2=np.zeros(ndays)
fieldrms_vec2=np.zeros(ndays)
ii=0
jj=0
for years in range(yearstart,yearend+1):
    print years
    resultdir1=resultbasedir+str(years)+'/'
    resultdir2=avisobasedir+str(years)+'/'
    filelist1=sorted(glob.glob(resultdir1+resultbasefile))
    filelist2=sorted(glob.glob(resultdir2+avisobasefile))
    
    for resfiles in filelist1:
        #print resfiles
        with netcdf.Dataset(resfiles) as nc:
            field=nc.variables['analyzed_field'][:]*100.
            fieldmean_vec[ii]=field.mean()
            fieldmin_vec[ii]=field.min()
            fieldmax_vec[ii]=field.max()
            ii+=1
    for resfiles in filelist2:
        with netcdf.Dataset(resfiles) as nc:
            field=nc.variables['sla'][:]*100.
            fieldmean_vec2[jj]=field.mean()
            fieldmin_vec2[jj]=field.min()
            fieldmax_vec2[jj]=field.max()
            jj+=1


# Create time to plot
time0=np.arange(0,ndays)
timeinit=calendar.timegm((yearstart,1,1,0,0,0))
time2plot=timeinit+time0*86400.

newxticks=np.zeros(nyears)
newxticklabels=[]
for yy,years in enumerate(range(yearstart,yearend+1)):
    newxticks[yy]=calendar.timegm((years,1,1,0,0,0))
    newxticklabels.append(str(years))

    
# Make some plots
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(time2plot+30.*86400,fieldmean_vec,label='DIVA')
plt.fill_between(time2plot+30.*86400, fieldmin_vec, fieldmax_vec, facecolor='blue', alpha=0.25)
plt.plot(time2plot,fieldmean_vec2,'r',label='AVISO')
plt.fill_between(time2plot, fieldmin_vec2, fieldmax_vec2, facecolor='red', alpha=0.25)
ax.set_xticks(newxticks[::2])
ax.set_xticklabels(newxticklabels[::2])
ax.set_xlim(calendar.timegm((yearstart,1,31,0,0,0)),calendar.timegm((yearend,12,31,0,0,0)))
ax.set_ylim(-40.,40.)
plt.legend(loc=2)
plt.xlabel('Time')
plt.title('SLA (cm)', fontsize=22)
plt.savefig(figdir+figname1, dpi=resolution, facecolor='w', edgecolor='w',
          transparent=True, bbox_inches='tight', pad_inches=0.)

ax.set_xlim(calendar.timegm((2006,1,1,0,0,0)),calendar.timegm((yearend,12,31,0,0,0)))
plt.savefig(figdir+figname2, dpi=resolution, facecolor='w', edgecolor='w',
          transparent=True, bbox_inches='tight', pad_inches=0.)
plt.close()


                              
    
