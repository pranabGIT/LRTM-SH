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

######################################
            # TO BE MODIFIED
######################################


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

# loading NCEP anom data 

cmappth = r"/media/pranab/Backup Plus/Backup_12_05_2021/NCEP_gh/"

flp = "hgtano250_NCEP_NCAR_1980JAN_2005DEC_daily_SH.nc"

flnamez = cmappth+"/"+flp

ncz = NetCDFFile(flnamez, 'r')
prtr = ncz.variables['HGT'][:,0,:,:]
tm1 = ncz.variables['TIME']
lon = ncz.variables['LONN71_73']
lat = ncz.variables['LAT']

#######################################################
tm2 = np.linspace(1,len(prtr),len(prtr))
###############################################################################
yr1, yr2 = 1980, 2005
yln = range(yr1,yr2)
yrind = 0
tsp = 90

[mzg,nzg,ozg] = np.shape(prtr)

pam = np.zeros((1,nzg,ozg))+99999.99
for iy in range(len(yln)): 
    yr = yln[iy]
    print ('years',yr, '-',str(yr+1))
    
    if (yr%4 == 0) or (yr%400 == 0):
        indst, indend = yrind+335, yrind+335+90 # 335 for until end of November, so that it starts from 1 Dec
        yrtot = 366
    else:
        indst, indend = yrind+334, yrind+334+90
        yrtot = 365
        
#     # FOR models with dy = 365
#     if dy == 365.0:
#         indst, indend = yrind+334, yrind+334+90
#         yrtot = 365

#     # FOR models with dy = 360
#     if dy == 360.0:
#         indst, indend = yrind+11*30, yrind+11*30+90 # 9*30=270 for until end of Sepember, 20 for October, so that it starts from 21 Oct ~ this gives us a lag of 40 days before DJF
#         yrtot = 360
    
    # extract the djf days daily data for #iy year
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
lon = np.array(lon)
lat = np.array(lat)

floutM = ('sig_mat_DJF_'+str(yr1)+'_'+str(yr2-1)+'_NCEP.mat')
print (floutM)


#np.savez(floutPzg, sig_mat=zg_sigmat, lon=lonz, lat=latz)
sio.savemat(floutM, {'sig_mat':sig_mat, 'lon':lon, 'lat':lat})
