#!/usr/bin/python

"""
# this code relates to our investigation of linear response 
in signals due to forcings in Maritime continent
(check emails from Adrian Mathews, dated: 14/12/2017)

#####################################
** (1) takes in daily cmap prec anom 
** (2) makes DJF prec anom composite during ENSO years (during 19800101-20051231 period)
** (3) calculates the projection of (2) on (1)
    ** computed the DJF frc_mat for (3)

** Dec centred, i.e., 2000-2009 starts with dec/2000-jan/2001-feb/2001


If the code doesn't run, or gives error 'Invalid Display Variable',
use the following line in the terminal before launching python:
$ export DISPLAY=localhost:0
"""

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


from IPython import get_ipython
ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')


from datetime import datetime
#######################################################

# loading TRMM anom data for 1998-2017

# DJF daily anomalies for all days/months/years - 1998-2017
# DJF daily anomalies for all days/months/years - 1998-2017
# cmappth = r"I:\CMAP-precip"
cmappth = r"/media/pranab/Backup Plus/Backup_12_05_2021/CMAP-precip"

flp = "CMAP_daily_anomaly_enh_01011980_31122005_TROPICS.nc"

flnamez = cmappth+"/"+flp

ncz = NetCDFFile(flnamez, 'r')
prtr = ncz.variables['PRECIP']
tm1 = ncz.variables['TIME']
lon = ncz.variables['LONN71_72']
lat = ncz.variables['LAT']


# flpr = cmappth + "/" + 'cmap_dailyanom_01011979_17122016_tropics.mat'

# prec = sio.loadmat(flpr)
# lat = prec["lat"]
# lon = prec["lon"]
# prtr = prec["prec"]
# tm1 = prec["tm"]
                             
tm2 = (tm1[:] - tm1[0]) + 1.0

[m,n,o] = np.shape(prtr)


#################################################################################
    # finds the DJF prec for all years: 1998-2017 #
yr1 = 1980
yr2 = 2005

yln = range(yr1, yr2)
yrind = 0
pryr = np.zeros((n, o, len(yln)))
for iy in range(len(yln)): 
    yr = yln[iy]
    print ('years',yr, '-',str(yr+1))
    if (yr%4 == 0) or (yr%400 == 0):
        indst, indend = yrind+335, yrind+335+90 # 335 for until end of November, so that it starts from 1 Dec
        # this gives us a lag of 40 days before DJF
        yrtot = 366
    else:
        indst, indend = yrind+334, yrind+334+90
        yrtot = 365

    # extract the djf (i.e., 90 days) daily data for #iy year
    pr = prtr[indst:indend,:,:]
    prav = np.nanmean(pr,0)
    #pdat = np.vstack((pdat,prav))

    tm3 = tm2[indst:indend]
    print ('time limits', tm3[0], tm3[-1])
    #print ('length of prec. DJF i.e. 90 days', len(pr))
 
    pryr[:,:,iy] = prav
    yrind = yrind + yrtot
    
####### ENSO COMPOSITE DJF #######
# extracts the composite prec for ENSO extreme years #

yln = np.array(yln)
 
# selected extreme years
# El Nino years
yr1 = [1982, 1987, 1991, 1997, 2002] # list for CMAP period
# La Nina events
# yr1 = [1988, 1995, 1998, 1999, 2000]; # list for CMAP period

pdat = np.zeros((n, o, 1))+99999.99
for iyr in range(len(yr1)):
    print(yr1[iyr])
    pp = np.where (yln == yr1[iyr])
    pp = pp[0]
    prc10 = pryr [:, :, pp]
    pdat = np.concatenate((pdat, prc10), axis=2)
    del prc10
    del pp

prc_cc = np.nanmean(pdat [:,:,1:], 2)

# floutP = 'CMAP_ElNino_compAnom_1980_2004DJF.npz'
# np.savez(floutP, prc_cc=prc_cc, lonc=lon, latc=lat)
     
###############################################################################
# compute the projection values and arrange as frc_mat

[m,n,o] = np.shape(prtr)


# find the projection of climate change on trmm - inner product
pfcc = []
for i in range(m):
        a = prtr[i,:,:]
        b = prc_cc # prec anom during El Nino
        prj = np.nansum(a*b)
        if np.isnan(prj) == 1:
            print ('inner product : ',prj)
        pfcc.append(prj)

asn = np.isnan(a)
N = n*o-np.sum(asn) # removing 'NaN' grid points near the boundary

# NaN values has to be removed, otherwise norm becomes NaN
b[np.isnan(b)] = 0 ## IMPORTANT: used to replace nan values with 0

fnrm = np.linalg.norm(b)

# normalising-1
pfcc = pfcc/(N*fnrm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
q75, q25 = np.percentile(pfcc, [75 ,25])
icc75 = np.where(pfcc>=q75)
icc25 = np.where(pfcc<=q25)
prtr75 = np.mean(prtr[icc75[0],:,:],0)
prtr25 = np.mean(prtr[icc25[0],:,:],0)
# plot prec ano corresponding to high/low/intermediate pfcc values
plt.close("all")
fig = plt.figure(figsize=(8,8))
[ln,lt] = np.meshgrid(lon,lat)
for ipl in range(3):
        if ipl == 0:
                ax1 = fig.add_subplot(311)
                pr = prc_cc
        if ipl == 1:
                ax1 = fig.add_subplot(312)
                pr = prtr75
        if ipl == 2:
                ax1 = fig.add_subplot(313)
                pr = prtr25

        parallels = np.arange(-20.,90.,20.)
        meridians = np.arange(0.,360.,90.)
        m = Basemap(llcrnrlon=0.,llcrnrlat=-20.,urcrnrlon=360.,urcrnrlat=20.,\
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='l',projection='merc',\
                    lat_0=0.,lon_0=180.,lat_ts=0.)
        clevsp = np.linspace(-5.,5.,11)
#        x,y = np.meshgrid(lon,lat)
        cs = m.contourf(ln,lt,pr,clevsp,cmap=plt.cm.RdBu,latlon=True,animated=True, extend = 'both')
        cbar = m.colorbar(cs,location='right',pad="10%")
        m.drawparallels(parallels,labels=[1,1,1,1],fontsize=7)
        m.drawmeridians(meridians,labels=[1,1,1,1],fontsize=7)
        m.drawcoastlines(color = 'gray', linewidth=0.5)
        
figname = 'pfccVdaily_ensoONcmap.png'

# plt.savefig(figname, dpi=400)


import sys
sys.exit()
##################################################################
# Extracting data for DJF+40 days lag
pfcc=scipy.signal.detrend(pfcc, axis=0)
#################################################################

print ('arranging prec data with lags as forcing matrix for finding correlations....')

yr1, yr2 = 1980, 2005
yln = range(yr1,yr2)

ilag = 40
tsp = 90

yrind = 0
pam = np.zeros((ilag+1))+99999.99
for iy in range(len(yln)): 
    yr = yln[iy]
    print ('years',yr, '-',str(yr+1))
    if (yr%4 == 0) or (yr%400 == 0):
        indst, indend = yrind+273+22, yrind+273+22+130 # 273 for until end of Sepember, 22 for October, so that it starts from 22 Oct
        # this gives us a lag of 40 days before DJF
        yrtot = 366
    else:
        indst, indend = yrind+272+22, yrind+272+22+130
        yrtot = 365

    # extract the djf+40 days daily data for #iy year
    pr = pfcc[indst:indend]

    tm3 = tm2[indst:indend]
    print ('time limits', tm3[0], tm3[-1])
    print ('length of prec. DJF+40 days', len(pr))
    # arrange prec data as forcing matrix for finding correlations
    for im in range(tsp):
           print (im)
           print (im+ilag)
           pa = pr[im:im+ilag+1]
           time2 = tm3[im:im+ilag+1]
           print (time2[0])
           print (time2[-1])
           print (len(time2))
           pam = np.vstack((pam,pa))

    yrind = yrind+yrtot

frc_mat = pam[1:,:] # this removes the first dummy row
print ('shape of frc_mat: ', np.shape(frc_mat))

floutM = ('Proj_cmapElNinoONcmap_frc_mat_lag'+str(ilag)+'_DJF_'+str(yr1)+'_'+str(yr2-1)+'.mat')


#print floutP
print (floutM)
#np.savez(floutP, frc_mat=frc_mat)
sio.savemat(floutM, {'frc_mat':frc_mat})
