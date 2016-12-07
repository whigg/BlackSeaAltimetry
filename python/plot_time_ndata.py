import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy  import *

filedir='/home/ctroupin/DIVA/BlackSea_Altimetry/diagnostics/'
filename='datanum_blacksea_daily.dat'
figdir='/home/ctroupin/DIVA/BlackSea_Altimetry/figures/'
figname=figdir+filename[0:-4]

# Load data from file
data=np.loadtxt(filedir+filename)
xmax=max(data[:,0])
ymax=max(data[:,1])
xmin=min(data[:,0])
ymin=min(data[:,1])


xx1=xmin+(xmax-xmin)/2.0
xx2=xmin+0.02*(xmax-xmin)
yy1=ymin+(ymax-ymin)/2.0
yy2=0.05*ymax

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(data[:,1]-data[:,1].min(),'ko-', ms=2, lw=1)


##ax.text(xx1,yy2,'Number of data points',ha='center',fontsize=18,fontname='Times New Roman')
##ax.text(xx2,yy1,'Computation time (s)',fontsize=18,fontname='Times New Roman')
##ax.set_xlim(xmin,xmax)
##ax.set_ylim(0,ymax)

#fig.autofmt_xdate()
#plt.savefig(figname, dpi=300, bbox_inches='tight')

plt.show()
plt.close()
