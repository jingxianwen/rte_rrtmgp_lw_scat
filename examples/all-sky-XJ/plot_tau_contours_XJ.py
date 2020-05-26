#=====================================================
#
#=====================================================
# os
import os

#import netCDF4
from netCDF4 import Dataset as netcdf_dataset

# cartopy
#import cartopy.crs as ccrs
#from cartopy.mpl.geoaxes import GeoAxes
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#from cartopy.util import add_cyclic_point

# matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors

# numpy
import numpy as np

# scipy
#from scipy import stats

# parameters
#from get_parameters import get_area_mean_min_max

fpath_ctl="./ref/rrtmgp-allsky-std-ice-ref-XJ.nc"         #control data
fpath_exp="./ref/rrtmgp-allsky-mc6-ice-ref-XJ.nc"  #experimental data

varnm="lw_tau_tot"      #the variable to be plotted
longname="Optical Depth"  

label_ctl="Standard"    #name of control data
label_exp="MC6_Scat"  #name of experimental data

file_ctl=netcdf_dataset(fpath_ctl,"r") 
file_exp=netcdf_dataset(fpath_exp,"r") 

dtctl=file_ctl.variables[varnm]
dtexp=file_exp.variables[varnm]

nbnd,ncol=dtctl.shape

var1=dtctl[:,:]
var2=dtexp[:,:]
vardiff=var2-var1
units="unit"

x=np.linspace(11.3,180,128)  #the Rei coordinate 
y=np.linspace(0.5,nbnd+0.5,nbnd)   #the band coordinate
#print(x)
#print(y)

#-------------------------------
# make the plot
#-------------------------------
fig=plt.figure(figsize=(6,8))
panel = [(0.15, 0.7, 0.7, 0.23),\
         (0.15, 0.4, 0.7, 0.23),\
         (0.15, 0.1, 0.7, 0.23),\
        ]  
ax1=fig.add_axes(panel[0])
ax2=fig.add_axes(panel[1])
ax3=fig.add_axes(panel[2])

levels1=np.linspace(1,10,10)
levels2=levels1
levels3=np.linspace(-4,4,9)

p1=ax1.contourf(x,y,var1,levels=levels1,cmap="bwr",extend="both")
p2=ax2.contourf(x,y,var2,levels=levels2,cmap="bwr",extend="both")
p3=ax3.contourf(x,y,vardiff,levels=levels3,cmap="bwr",extend="both")

cbax = fig.add_axes((panel[0][0] + 0.72, panel[0][1] + 0.0354, 0.0326, 0.1792))
cbar = fig.colorbar(p1, cax=cbax, ticks=levels1)
cbax = fig.add_axes((panel[1][0] + 0.72, panel[1][1] + 0.0354, 0.0326, 0.1792))
cbar = fig.colorbar(p2, cax=cbax, ticks=levels2)
cbax = fig.add_axes((panel[2][0] + 0.72, panel[2][1] + 0.0354, 0.0326, 0.1792))
cbar = fig.colorbar(p3, cax=cbax, ticks=levels3)

ax1.set_title(longname+" ("+label_ctl+")",fontsize=12)
ax2.set_title(longname+" ("+label_exp+")",fontsize=12)
ax3.set_title(longname+" ("+label_exp+"-"+label_ctl+")",fontsize=12)
ax1.set_ylabel("Bands",fontsize=12)
ax2.set_ylabel("Bands",fontsize=12)
ax3.set_ylabel("Bands",fontsize=12)
#ax1.set_xlabel("Rei (micron)",fontsize=14)
#ax2.set_xlabel("Rei (micron)",fontsize=14)
ax3.set_xlabel("Rei (micron)",fontsize=12)

#plt.savefig(figure_name+".png",dpi=(200))
plt.show()

exit()
