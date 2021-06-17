"""
# this code relates to our investigation of linear response 
in signals due to forcings in Maritime continent
(check emails from Adrian Mathews, dated: 14/12/2017)

#####################################
Script to compute signal matrix of daily ZG/Prec for DJF
over North Pacific :: signal

FOR both HISTORICAL (1985-2005) and RCP 8.5 scenarios (2079-2099)

(input: uses already detrended but unfiltered anomaly data for all months
created by 'anomaly***.py' code)

[*****modify for other seasons*****]

** Finally, producing 'm x n x o' matrix with m -> time dim, n -> lat dim, o -> lon dim

m = 90 days of DJF x no.s of years corresponding to the decade

*************************************
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

cmippth = r"/media/pranab/Backup Plus/Backup_12_05_2021/LRTM-SH-PD-ProcessedData"

# CCSM4
#flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_CCSM4_amip_19800101_20051231.npz'; modl = 'CCSM4_amip'; dy = 365.0; print ('Model chosen :: CCSM4')

#HadGEM2A
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_HadGEM2A_amip_19800101_20051230.npz'; modl = 'HadGEM2A_amip'; dy = 360.0; print ('Model chosen :: HadGEM2A') 

#MIROC5
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_MIROC5_amip_19800101_20051231.npz'; modl = 'MIROC5_amip'; dy = 365.0; print ('Model chosen :: MIROC5')

#IPSLcm5aMR
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_IPSLcm5aMR_amip_19800101_20051231.npz'; modl = 'IPSLcm5aMR_amip'; dy = 365.0; print ('Model chosen :: IPSLcm5aMR')

# NorESM1
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_NorESM1_amip_19800101_20051231.npz'; modl = 'NorESM1_amip'; dy = 365.0; print ('Model chosen :: NorESM1')

# GFDL-CM3-amip: USES GREGORIAN CALENDAR
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_GFDLCM3_amip_19800101_20051231.npz'; modl = 'GFDLCM3_amip'; dy = 365.0; print ('Model chosen :: GFDL-CM3');

# MPI-ESM-MR-amip: USES GREGORIAN CALENDAR
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_MPIesmMR_amip_19800101_20051231.npz'; modl = 'MPIesmMR_amip'; dy = 365.0; print ('Model chosen :: MPI-ESM-MR');

# ACCESS1-3-amip: USES GREGORIAN CALENDAR
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_ACCESS1-3_amip_19800101_20051231.npz'; modl = 'ACCESS1-3_amip'; dy = 365.0; print ('Model chosen :: ACCESS1-3');

# MRI-CGCM3-amip: USES GREGORIAN CALENDAR
# flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_MRI-CGCM3_amip_19800101_20051231.npz'; modl = 'MRI-CGCM3_amip'; dy = 365.0; print ('Model chosen :: MRI-CGCM3');


# NorESM1-amip: USES GREGORIAN CALENDAR
flpr = cmippth + "/" + 'zg250_ano_daily_detrend_NoFilt_SH_NorESM1_amip_19800101_20051231.npz'; modl = 'NorESM1_amip'; dy = 365.0; print ('Model chosen :: NorESM1');



prfile = np.load(flpr)

prtr = prfile['zgano']
[mnz,nzg,ozg] = np.shape(prtr)
lon = prfile['lon']
lat = prfile['lat']
[lnx, lty] = np.meshgrid(lon, lat)
#######################################################

tm2 = np.linspace(1,len(prtr),len(prtr))
###############################################################################
yr1, yr2 = 1980, 2005
yln = range(yr1,yr2)
yrind = 0
tsp = 90


pam = np.zeros((1,nzg,ozg))+99999.99
for iy in range(len(yln)): 
    yr = yln[iy]
    print ('years',yr, '-',str(yr+1))
    
    #if (yr%4 == 0) or (yr%400 == 0):
        #indst, indend = yrind+273+22, yrind+273+22+130 # 273 for until end of Sepember, 22 for October, so that it starts from 22 Oct
        # this gives us a lag of 40 days before DJF
        #yrtot = 366
    #else:
    
    # FOR models with dy = 365
    if dy == 365.0:
        indst, indend = yrind+334, yrind+334+90
        yrtot = 365

    # FOR models with dy = 360
    if dy == 360.0:
        indst, indend = yrind+11*30, yrind+11*30+90 # 9*30=270 for until end of Sepember, 20 for October, so that it starts from 21 Oct ~ this gives us a lag of 40 days before DJF
        yrtot = 360
    
    # extract the djf+40 days daily data for #iy year
    pr = prtr[indst:indend,:,:]
    tm3 = tm2[indst:indend]
    
    if len(tm3) == 90:
        pam = np.vstack((pam,pr))
    else:
        print ('CHECK IF IT IS AN ERROR OR IF IT IS SIMPLY THE LAST YEAR (WHICH WONT HAVE 90 days in DJF)')
        # this gives us a lag of 40 days before DJF

    
    print ('time limits', tm3[0], tm3[-1])
    print ('length of zg DJF (90) days', len(pr))
    
    yrind = yrind+yrtot

sig_mat = pam[1:,:,:] # this removes the first dummy row

print ('shape of sig_mat: ', np.shape(sig_mat))

# import sys
# sys.exit()

floutM = ('sig_mat_DJF_'+str(yr1)+'_'+str(yr2-1)+'_'+modl+'.mat')
print (floutM)

#
#np.savez(floutPzg, sig_mat=zg_sigmat, lon=lonz, lat=latz)
sio.savemat(floutM, {'sig_mat':sig_mat, 'lon':lon, 'lat':lat})
