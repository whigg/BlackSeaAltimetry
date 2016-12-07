#!/usr/bin/env python

# diva_plot_mesh.py
#
# Plot the finite-element mesh
#
# http://modb.oce.ulg.ac.be/mediawiki/index.php/Diva_python
#------------------------------------------------------------------------------

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
from matplotlib.path import Path
import matplotlib.patches as patches

# Clean 
os.system('clear')

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'l'

# File and directory names
meshdir='/home/ctroupin/DIVA/BlackSea_Altimetry/mesh/L1/'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/mesh/'
meshfile = 'mesh.dat'
meshtopofile = 'meshtopo.dat'
coastdir='/home/ctroupin/DataOceano/Coastlines/'
coastfile='blacksea_coast2.dat'
datadir='/home/ctroupin/DataOceano/AVISO/BlackSea/Unfiltered/2-mission_concat/'
datafile='dt_blacksea_c2ennj2_sla_vxxc_20120328.dat'

cmap=plt.cm.RdYlBu_r

# Limits of the colorbar
vmin,vmax,dvar = -20.0,20.0,5.0
bounds2 = np.arange(vmin,vmax+0.0001,dvar)
norm = colors.Normalize(vmin=vmin,vmax=vmax)



# Figure name
figbasename = 'blacksea_mesh_data'
figtype='.png'
figname=figdir+figbasename+figtype

# Region of interest        
lonmin=27
lonmax=42
latmin=40
latmax=48

# Spacing between x/y labels
dlon = 3
dlat = 2

# Unit of the variable to plot
unitname = ' '

#------------------------------------------------------------------------------------

# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

# ----------------------------------------------------------------------------

    
# Load mesh information
datamesh = np.loadtxt(meshdir+meshtopofile)
nnodes = int(datamesh[0])
ninterfaces = int(datamesh[1])
nelements = int(datamesh[2])
ntotal = int(datamesh.sum())

# Load mesh nodes
meshnodes = np.genfromtxt(meshdir+meshfile,skip_footer=nelements+ninterfaces)
meshnodes = np.fromstring(meshnodes)

# Load mesh elements
meshelements = np.genfromtxt(meshdir+meshfile,skip_header=nnodes+ninterfaces)
meshelements = np.fromstring(meshelements)
meshelements = np.int_(meshelements)

# Extract node coordinates
xnode=meshnodes[np.arange(1,nnodes*3,3)]
ynode=meshnodes[np.arange(2,nnodes*3,3)]

# Indices of the elements
i=np.transpose(range(0,nnodes));
i1=meshelements[np.arange(0,nelements*6,6)]-1
i2=meshelements[np.arange(2,nelements*6,6)]-1
i3=meshelements[np.arange(4,nelements*6,6)]-1

# Make the plot
fig=plt.figure()
ax = fig.add_subplot(111)

m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                lat_ts=0.5*(lonmin+lonmax),\
                resolution=basemap_resolution)
m.ax=ax

# Load extracted coastline
valexc=-999.
lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
lonc=np.ma.masked_where(lonc==valexc,lonc)
latc=np.ma.masked_where(latc==valexc,latc)
xc,yc=m(lonc,latc)

# Project nodes
xnode,ynode=m(xnode,ynode)
xnode[xnode==1e+30]=np.nan
ynode[ynode==1e+30]=np.nan

# Loop on the elements and patch
# (can probably be optimized if the loop is removed)

for j in range(0,nelements):
    verts = [(xnode[i1[j]],ynode[i1[j]]),\
             (xnode[i2[j]],ynode[i2[j]]),\
             (xnode[i3[j]],ynode[i3[j]]),\
             (xnode[i1[j]],ynode[i1[j]]) ]
    path = Path(verts)
    patch = patches.PathPatch(path, facecolor='none', lw=1)
    m.ax.add_patch(patch)

# load data
lon,lat,field,=np.loadtxt(datadir+datafile,usecols=(0,1,2),unpack='True')
field*=100.
x,y = m(lon, lat)
        
# Set axis limits 
m.ax.set_xlim(lonmin,lonmax)
m.ax.set_ylim(latmin,latmax)

# Add grid, coastline and continent
m.plot(xc,yc,'k',lw=0.5)
#m.drawcoastlines(ax=ax)
#m.fillcontinents(color='black', ax=ax)
m.bluemarble(ax=ax,alpha=0.85)
scat = m.scatter(x, y, c=field,s=15,zorder=2, edgecolor='none',norm=norm,cmap=cmap)

meridians=np.arange(lonmin,lonmax,dlon)
parallels=np.arange(latmin,latmax,dlat)
m.drawmeridians(meridians,labels=[0, 0, 0, 1],fontsize=18,fontname='Times New Roman')
m.drawparallels(parallels,labels=[1, 0, 0, 0],fontsize=18,fontname='Times New Roman')

plt.title('2012-03-28',fontsize=22,fontname='Times New Roman')
# Export figure and display it
plt.savefig(figname, dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
plt.close()

