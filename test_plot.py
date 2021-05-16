#!/usr/bin/python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import glob, os
from netCDF4 import Dataset as NetCDFFile
import math
import os
os.environ['PROJ_LIB'] = r"C:\Users\Pranab\Anaconda3\Library\share"
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap, cm
matplotlib.pyplot.switch_backend('agg')
from matplotlib.patches import Ellipse, Polygon
import scipy.io as sio
import scipy.signal
from datetime import datetime
from scipy.stats.stats import pearsonr
import scipy.interpolate
from IPython import get_ipython
ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')
from datetime import datetime
#######################################################
# DATA
cmippth = r"I:\LRTM-SH-PD-ProcessedData"
flpr = cmippth + "\\" + 'pr_ano_daily_detrend_NoFilt_NH_CCSM4_amip_19800101_20051231.npz'
prfile = np.load(flpr)
prtr = prfile['prano']
lon = prfile['lon']
lat = prfile['lat']
[lnx, lty] = np.meshgrid(lon, lat)
#######################################################
# loads the composite prec ano for ENSO extreme years #
ifl = 1
flc = 'CMAP_ElNino_compAnom.npz'; flo = 'ElNino'
flEN = np.load(flc)
prEN = flEN['prc_cc']
lonc = flEN['lonc']
latc = flEN['latc']

#y,x=np.meshgrid(lonc,latc)
y,x=np.meshgrid(latc,lonc)

# prtr --> prEN
# flattened lon/lat for regridding
#points = np.array((lnx.flatten(),lty.flatten())).T
#prc_cc = scipy.interpolate.griddata(points, prtr[0,:,:].flatten(), (x,y),method='cubic')

# prEN --> prtr
# flattened lon/lat for regridding
points = np.array((x.flatten(),y.flatten())).T

#prtr = griddata(points, prano, (lnx,lty))

#prc_cc=scipy.interpolate.griddata((x.flatten(),y.flatten()), prEN.flatten(), (lnx,lty),method='cubic')
prc_cc = scipy.interpolate.griddata(points, prEN.flatten() , (lnx,lty),method='cubic')
#######################################################
# PLOTS
parallels = np.arange(-20.,90.,20.)
meridians = np.arange(0.,360.,90.)
plt.close("all")
fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(211)
m = Basemap(llcrnrlon=0.,llcrnrlat=-20.,urcrnrlon=360.,urcrnrlat=20.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=0.,lon_0=180.,lat_ts=0.)

clevsp = np.linspace(-5.,5.,11)
#        x,y = np.meshgrid(lon,lat)
cs = m.contourf(x,y,prEN,clevsp,cmap=plt.cm.RdBu,latlon=True,animated=True, extend = 'both')

cbar = m.colorbar(cs,location='right',pad="10%")

m.drawparallels(parallels,labels=[1,1,1,1],fontsize=7)
m.drawmeridians(meridians,labels=[1,1,1,1],fontsize=7)
m.drawcoastlines(color = 'gray', linewidth=0.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ax2 = fig.add_subplot(212)
m = Basemap(llcrnrlon=0.,llcrnrlat=-20.,urcrnrlon=360.,urcrnrlat=20.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=0.,lon_0=180.,lat_ts=0.)

clevsp = np.linspace(-5.,5.,11)
#        x,y = np.meshgrid(lon,lat)
cs = m.contourf(lnx,lty,prc_cc,clevsp,cmap=plt.cm.RdBu,latlon=True,animated=True, extend = 'both')
cbar = m.colorbar(cs,location='right',pad="10%")

m.drawparallels(parallels,labels=[1,1,1,1],fontsize=7)
m.drawmeridians(meridians,labels=[1,1,1,1],fontsize=7)
m.drawcoastlines(color = 'gray', linewidth=0.5)

