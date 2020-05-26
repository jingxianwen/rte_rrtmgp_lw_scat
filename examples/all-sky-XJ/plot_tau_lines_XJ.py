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
ibnd= min(3,nbnd-1)     #which band to plot? 

var1=dtctl[ibnd,:]
var2=dtexp[ibnd,:]
longname=longname+" (band "+str(ibnd+1)+")"
units="unit"

x=np.linspace(11.3,180,128)  #the Rei coordinate 
#print(x)

#-------------------------------
# make the plot
#-------------------------------
fig=plt.figure(figsize=(8,6))
ax1=fig.add_axes([0.12,0.24,0.8,0.6])

ax1.plot(x,var1,color="k",lw=2,ls="-",label=label_ctl)
ax1.plot(x,var2,color="r",lw=2,ls="-",label=label_exp)
ax1.legend(fontsize=12)
ax1.set_title(longname,fontsize=14)
ax1.set_ylabel(units,fontsize=14)
ax1.set_xlabel("Rei (micron)",fontsize=14)
ax1.set_xlim(min(x),max(x))

#plt.savefig(figure_name+".png",dpi=(200))
plt.show()

exit()
