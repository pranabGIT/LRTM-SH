#!/usr/bin/python

"""
# this code relates to our investigation of linear response 
in signals due to forcings in Maritime continent
(check emails from Adrian Mathews, dated: 14/12/2017)

#####################################
** (1) loads daily prec anom from various CMIP models to give composite DJF for El Nino years
** (2) loads the anom of daily prec from CMAP
** I M P O R T A N T: regrids (1) on (2) grid
** (3) calculates the projection of (1) on (2)
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
import scipy.interpolate


from IPython import get_ipython
ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')


from datetime import datetime
#######################################################

# loading CMIP anom data 

#cmippth = '/media/pranab/STORAGE5/LRTM-SH-PD-ProcessedData'
#cmippth = '/media/pranab/Backup Plus'
cmippth = '/media/pranab/Backup Plus/Backup_12_05_2021/LRTM-SH-PD-ProcessedData'

# CCSM4
flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_CCSM4_amip_19800101_20051231.npz'; modl = 'CCSM4_amip'; dy = 365.0; print ('Model chosen :: CCSM4')

#HadGEM2A
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_HadGEM2A_amip_19800101_20051230.npz'; modl = 'HadGEM2A_amip'; dy = 360.0; print ('Model chosen :: HadGEM2A') 

#MIROC5
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_MIROC5_amip_19800101_20051231.npz'; modl = 'MIROC5_amip'; dy = 365.0; print ('Model chosen :: MIROC5')

#IPSLcm5aMR
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_IPSLcm5aMR_amip_19800101_20051231.npz'; modl = 'IPSLcm5aMR_amip'; dy = 365.0; print ('Model chosen :: IPSLcm5aMR')

#GFDL-CM3
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_GFDLCM3_amip_19800101_20051231.npz'; modl = 'GFDLCM3_amip'; dy = 365.0; print ('Model chosen :: GFDLCM3')

#MPI-ESM-MR
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_MPIesmMR_amip_19800101_20051231.npz'; modl = 'MPIesmMR_amip'; dy = 365.0; print ('Model chosen :: MPI-ESM-MR')

#ACCESS1-3
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_ACCESS1-3_amip_19800101_20051231.npz'; modl = 'ACCESS1-3_amip'; dy = 365.0; print ('Model chosen :: ACCESS1-3')

#MRI-CGCM3
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_MRI-CGCM3_amip_19800101_20051231.npz'; modl = 'MRI-CGCM3_amip'; dy = 365.0; print ('Model chosen :: MRI-CGCM3')

#NorESM1
# flpr = cmippth + "/" + 'pr_ano_daily_detrend_NoFilt_NH_NorESM1_amip_19800101_20051231.npz'; modl = 'NorESM1_amip'; dy = 365.0; print ('Model chosen :: NorESM1')

prfile = np.load(flpr)

prtr = prfile['prano']
lon = prfile['lon']
lat = prfile['lat']
tm2 = np.linspace(1,len(prtr),len(prtr))
[m,n,o] = np.shape(prtr)
# y,x=np.meshgrid(lat,lon)
x,y=np.meshgrid(lon,lat)
points = np.array((x.flatten(),y.flatten())).T

######### FIND THE MODEL COMPOSITE ENSO PREC #########
# Step-1: finds the DJF prec for all years: 1998-2017 #

yr1, yr2 = 1980, 2005
yln = range(yr1, yr2)
yln = np.array(yln)

yrind = 0
pryr = np.zeros((len(yln), n, o))
for iy in range(len(yln)): 
    yr = yln[iy]

    print ('year :: ', yr)
        # FOR models with dy = 365
    if dy == 365.0:
        indst, indend = yrind+334, yrind+334+90
        yrtot = 365

    # FOR models with dy = 360
    if dy == 360.0:
        indst, indend = yrind+11*30, yrind+11*30+90 # 9*30=270 for until end of Sepember, 20 for October, so that it starts from 21 Oct ~ this gives us a lag of 40 days before DJF
        yrtot = 360
   
    # extract the djf (i.e., 90 days) daily data for #iy year
    pr = prtr[indst:indend,:,:]
    prav = np.nanmean(pr,0)
    #pdat = np.vstack((pdat,prav))

    tm3 = tm2[indst:indend]
    print ('time limits', tm3[0], tm3[-1])
    #print ('length of prec. DJF i.e. 90 days', len(pr))
 
    pryr[iy,:,:] = prav
    yrind = yrind + yrtot
    
# Step-2: extracts the composite DJF prec for ENSO extreme years #
# El Nino years
# yr1 = [1982, 1987, 1991, 1997, 2002]; flo = 'ElNino'; # list for CMAP period

# La Nina events
yr1 = [1988, 1995, 1998, 1999, 2000]; flo = 'LaNina'; # list for CMAP period

pdat = np.zeros((1, n, o))+99999.99
for iyr in range(len(yr1)):
    
    pp = np.where (yln == yr1[iyr])
    pp = pp[0]
    prc10 = pryr [pp, :, :]
    pdat = np.concatenate((pdat, prc10), axis=0)
    del prc10
    del pp

prc_cc = np.nanmean(pdat [1:,:,:], 0)

#######################################################
# loading CMAP/Obs daily anomalies for all days/months/years - 1980-2005
# cmappth = r"I:\CMAP-precip"
cmappth = r"/media/pranab/Backup Plus/Backup_12_05_2021/CMAP-precip"
flp = "CMAP_daily_anomaly_enh_01011980_31122005_TROPICS.nc"
flnamez = cmappth+"/"+flp

ncz = NetCDFFile(flnamez, 'r')
prtr = ncz.variables['PRECIP']
tm1 = ncz.variables['TIME']
lonc = ncz.variables['LONN71_72']
lonc = np.array(lonc)
latc = ncz.variables['LAT']
latc = np.array(latc)
tm2 = (tm1[:] - tm1[0]) + 1.0
[m,n,o] = np.shape(prtr)

########################################################
# compute the projection values and arrange as frc_mat

# # flattened lon/lat for regridding regridding later (cmap grid to model grid)
[lnx, lty] = np.meshgrid(lonc, latc)
prc_cc = scipy.interpolate.griddata(points, prc_cc.flatten(), (lnx,lty), method='cubic')

[m,n,o] = np.shape(prtr)


# # find the projection of model enso comp on cmap - inner product
pfcc = []
for i in range(m):
         a = prtr[i,:,:]
#          a = scipy.interpolate.griddata(points, prtr[i,:,:].flatten(), (lnx,lty),method='cubic')
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

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

q75, q25 = np.percentile(pfcc, [75 ,25])

icc75 = np.where(pfcc>=q75)
icc25 = np.where(pfcc<=q25)

prtr75 = np.mean(prtr[icc75[0],:,:],0)
prtr25 = np.mean(prtr[icc25[0],:,:],0)

# plot prec ano corresponding to high/low/intermediate pfcc values
plt.close("all")
fig = plt.figure(figsize=(8,8))
#[ln,lt] = np.meshgrid(lon,lat)

for ipl in range(3):
        if ipl == 0:
                ax1 = fig.add_subplot(311)
                pr = prc_cc
                plt.title(modl)
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
        cs = m.contourf(lnx,lty,pr,clevsp,cmap=plt.cm.RdBu,latlon=True,animated=True, extend = 'both')
        cbar = m.colorbar(cs,location='right',pad="10%")

        m.drawparallels(parallels,labels=[1,1,1,1],fontsize=7)
        m.drawmeridians(meridians,labels=[1,1,1,1],fontsize=7)
        m.drawcoastlines(color = 'gray', linewidth=0.5)

figname = 'pfccVdaily_modlONcmap'+modl[0:4]+'.png'
plt.savefig(figname, dpi=400)

##################################################################
# Extracting data for DJF+40 days lag
pfcc=scipy.signal.detrend(pfcc, axis=0)
#################################################################

print ('arranging prec data with lags as forcing matrix for finding correlations....')

tm2 = np.linspace(1,len(prtr),len(prtr))
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

floutM = ('Proj_modlONcmap'+flo+'_frc_mat_lag'+str(ilag)+'_DJF_'+str(yr1)+'_'+str(yr2-1)+'_'+modl+'.mat')


#print floutP
print (floutM)


#import sys
#sys.exit()

#np.savez(floutP, frc_mat=frc_mat)
sio.savemat(floutM, {'frc_mat':frc_mat})

