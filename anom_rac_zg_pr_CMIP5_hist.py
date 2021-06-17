"""
- code finds a common period for prec and zg data (for hist period)
prepares the anomaly by removing annual cycle (smoothed with fft)
for that period, these wil be used to extract sig_mat and frc_mat

- saves the data for further use

- can be used for all cmip5 data on JASMIN

- Linear trend NOT removed from anomaly

- IMPORTANT:: CHECK FOR PI CONTROL OR FUTURE RUNS!!!!!
Pranab, Aug-2018

"""
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
###################################################
# CHOOSE THE MODEL #
###################################################
print ('CHOOSE MODEL ----')
print ('IMPORTANT:: CHECK FOR HIST OR PI CONTROL OR FUTURE RUNS!!!!!')

###################################################
        ########### YEAR RANGE ###########
yr1 = 1980; yr2 = 2005;

###################################################
# HIST
###################################################

# PROBLEMATIC and messed up Models (files do not have continuous time axis)

# BCC-CSM1-1 (files do not have continuous time axis)
#modnm = '/media/pranab/STORAGE5/cmip5-daily/BCC-CSM1-amip-hist-daily'; modl = 'BCC-CSM1_amip'; dy = 365.0; yrst = 1979.0; prc = 'pr'; time_offset = 0; print ('Model chosen :: BCC-CSM1'); g = 1;


# CMCC-CM (files do not have continuous time axis)
#modnm = "/media/pranab/STORAGE5/cmip5-daily/CMCC-CM-amip-hist-daily"; modl = 'CMCC-CM_amip'; dy = 365.0; yrst = 1979.0; prc = 'pr'; time_offset = 0; print ('Model chosen :: CMCC-CM'); g = 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOLLOWING MODELS ARE FINE

# CCSM4
modnm = "/media/pranab/Backup Plus/Backup_12_05_2021/cmip5-daily/CCSM4-amip-hist-daily"; modl = 'CCSM4_amip'; dy = 365.0; yrst = 1979.0; prc = 'prc'; time_offset = 0; print ('Model chosen :: CCSM4'); g = 0;

# MIROC5
#modnm = "/media/pranab/Backup Plus/Backup_12_05_2021/cmip5-daily/MIROC5-amip-hist-daily"; modl = 'MIROC5_amip'; dy = 365.0; yrst = 1979.0; prc = 'prc'; time_offset = 0; print ('Model chosen :: MIROC5'); g = 0;


# HadGEM2-ES :: Since original timeseries starts from 19780901, time_offset of 120 (4*30) is to be subtracted so that the time variable starts from 0.5 (i.e., from 19790101)
#modnm = "/media/pranab/Backup Plus/Backup_12_05_2021/cmip5-daily/HadGEM2-A-amip-hist-daily"; modl = 'HadGEM2A_amip'; dy = 360.0; yrst = 1979.0; prc = 'prc'; time_offset = 120; print ('Model chosen :: HadGEM2-A'); g = 0;


# IPSL-CM5A-MR

#modnm = "/media/pranab/Backup Plus/Backup_12_05_2021/cmip5-daily/IPSL-CM5A-MR-amip-hist-daily"; modl = 'IPSLcm5aMR_amip'; dy = 365.0; yrst = 1950.0; prc = 'prc'; time_offset = 0; print ('Model chosen :: IPSLcm5aMR'); g = 0;


# NorESM1
#modnm = '/media/pranab/Backup Plus/Backup_12_05_2021/cmip5-daily/NorESM1-M-amip-hist-daily'; modl = 'NorESM1_amip'; dy = 365.0; yrst = 1979.0; prc = 'pr'; time_offset = 0; print ('Model chosen :: NorESM1'); g = 0;


# Gregorian

# GFDL-CM3-amip: USES GREGORIAN CALENDAR
#modnm = "/media/pranab/STORAGE5/cmip5-daily/GFDL-CM3-amip-hist-daily"; modl = 'GFDLCM3_amip'; dy = 365.0; yrst = 1979.0; prc = 'prc'; time_offset = 0; print ('Model chosen :: GFDL-CM3'); g = 1;


# MPI ESM MR :: Gregrian
#modnm = "/media/pranab/STORAGE5/cmip5-daily/MPI-esm_MR-amip-hist-daily"; modl = 'MPIesmMR_amip'; dy = 365.0; yrst = 1979.0; prc = 'prc'; time_offset = 0; print ('Model chosen :: MPI ESM MR'); g = 1;


# MRI-CGCM3 :: Gregorian
#modnm = '/media/pranab/STORAGE5/cmip5-daily/MRI-CGCM3-amip-hist-daily'; modl = 'MRI-CGCM3_amip'; dy = 365.0; yrst = 1979.0; prc = 'pr'; time_offset = 0; print ('Model chosen :: MRI-CGCM3'); g = 1;


# ACCESS1-3 [yrst=1979 means time index for 19790101 is 0.5] :: Gregorian
#modnm = '/media/pranab/STORAGE5/cmip5-daily/ACCESS1-3-amip-hist-daily'; modl = 'ACCESS1-3_amip'; dy = 365.0; yrst = 1978.0; prc = 'pr'; time_offset = 722084; print ('Model chosen :: ACCESS1-3'); g = 1;


############### SEGMENT - 1 ################
# find common starting time step

    ######### zg #########
cpathz = modnm + "/" + "zg"

os.chdir(cpathz)

print ('###### list zg files ######')
#flz = []
#for file in glob.glob("*.nc"):
#    flz.append(file)
#    print(file)
    
flz = sorted(os.listdir(cpathz))

os.chdir('/home/pranab/Documents/LRTM-SH-PD')

flnamez = cpathz+"/"+flz[0]
ncz = NetCDFFile(flnamez, 'r')
tmz = ncz.variables['time'][0]
print ('1st zg file : ', flz[0])
print ('1st zg time step : ', tmz)

      ########## prec ###########

cpathp = modnm + "/" + "prc"

os.chdir(cpathp)

print ('###### list prec files #######')

flp = sorted(os.listdir(cpathp))

os.chdir('/home/pranab/Documents/LRTM-SH-PD')

flnamep = cpathp+"/"+flp[0]
ncp = NetCDFFile(flnamep, 'r')
tmp = ncp.variables['time'][0]
print ('1st pr file : ', flp[0])
print ('1st pr time step : ', tmp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find the common starting time step
# if prc starts at 1-jan 230 and zg starts at 1-jan, 430 ...
# take the starting point of data as 1-jan, 430
t11 =  max(tmp, tmz)
print ('common starting time step :: ', t11)

if tmp == tmz:
    print ('both zg and pr start from same timestep')
else:
    print ('zg and pr start from different time step:: S T R A N G E ! ! ! ! !')

tmp = tmp - time_offset;
tmz = tmz - time_offset;

t11 =  max(tmp, tmz)
print ('modified/offset starting time step :: ', t11)


# First and Last time step for whole period
t1 = (yr1-yrst)*dy+0.5 # Jan 1, 1980
print ('1st time step of 25-yr period (1980-2005):: ', t1)
t2 = t1+ (yr2-yr1+1)*dy -1 # Dec 31/30, 2005

# fllwoing step includes all the 29 Febs for Leap years (during 1980-2005) in models with gregorian calendar
if g == 1:
    t2 = t2+7
print ('last time step of 25-yr period (1980-2005):: ', t2)

print ('year1 of 25-yr period :: ', yr1)
print ('dec-31 of year25 of 25-yr period :: ', yr2)

print ('#################### S T A R T ####################')

print ('#################### EXTRACT 1980-2005 DATA from ORIGINAL FILES ####################')
# strating year
#yr1 = (t11-0.5)/dy+yrst
#print ('starting year :: ', yr1)

# find the last time step corresponding to 2005-12-31
#t2 = t11 + (2005-yr1)*dy -1  # upto last day of Dec, 2005
#t1 = t2 - 25*dy + 1 # 1 is added to get to jan -1

#print ('jan-1 of 1st year of the dataset :: ', t11)
#print ('starting timestep of last-20 year (i.e., jan-1 of year-1 of extracted 25 year period - 1980-2005) :: ', t1)
#print ('dec-31 of last year (i.e., of 2005) :: ', t2)
#print ('no. of years covered : ', (t2-t1+1)/dy)
#yr1_20 = (t1-0.5+1)/dy+yrst
#print ('year1 of 25-yr period :: ', yr1_20)
#yr2_20 = (t2-0.5+1)/dy+yrst # need to go to dec-1 (by adding 1) to calculate the year
#print ('dec-31 of year25 of 25-yr period :: ', yr2_20)

#########
### filename: should be manually put!!!!!!!!!!!!!!!!!!!******************
# storing the first & last date of the common period of data 
fyr1 = flp[0][-16:-12]
fyr2 = flp[0][-7:-3]

fyr = str(int(yr1))+fyr1+'_'+str(int(yr2))+fyr2

#print ('year extension for OP file: ', fyr)


################### SEGMENT - 2 #####################
#### extract pr & zg data for common time period ####
#####################################################

#### precipitation time slice ####
print ('#### extract PRECIPITATION #####')

# extract lat range and dim
flnamep = cpathp+"/"+flp[0]
ncp = NetCDFFile(flnamep, 'r')
lat1 = ncp.variables['lat'][:]
lonp = ncp.variables['lon'][:]
ilt = np.where ((lat1>=-20) & (lat1<=20))
latp = lat1[ilt[0]]
prtest = ncp.variables[prc][0,ilt[0],:]
[npr,opr] = np.shape(prtest)
del prtest

pr = np.zeros((1,npr,opr))
gdel = []; g1 = 0; g2 = 0;

for ip in range(len(flp)): # loop over all the prec files
    flpr = cpathp+"/"+flp[ip]
    print (flp[ip])
    nc = NetCDFFile(flpr, 'r')
    tim = nc.variables['time']
    tim = tim[:] - time_offset
    tim1 = tim[0]
    tim2 = tim[-1]
    print ('1st & last time steps in the file:: ', tim1, tim2)
    
    # find the time slices we need from each file

    if tim2 < t1:
        continue # go to next iteration
    elif tim1<=t1:
        int1 = np.where(tim[:] == t1)
        int1 = int1[0][0] # tuple to integer
    elif tim1>t1:
        int1 = 0

    if tim2 <= t2:  # t2 == 30 Nov of the final year (2045 here)
        int2 = len(tim)-1
    elif tim2 > t2-1:
        int2 = np.where(tim[:] == t2)
        int2 = int2[0][0] # tuple to integer
    
    print ('int1 and int2 - indices starting from Jan1-1980 :: ', int1,int2)
    print ('first and last time steps, starting from Jan1-1980 :: ', tim[int1], tim[int2])
# multiply by 86400 to convert from (kg m-2 s-1)  to (mm/day)
    pr1 =  nc.variables[prc][int1:int2+1, ilt[0], :]*86400 # timeslice up to (t2-1)//last index is intm in actuality i.e., a[0:5] will give a[0], a[1], a[2], a[3], a[4] i.e., 5 elements
    print ('size of extracted data from this file ::', np.shape(pr1))
    
        
    pr = np.vstack((pr,pr1))
    del pr1

    if tim2 >= t2-1:
        break
       
    
pr = pr[1:,:,:]
print ('size of extracted prec file (should be same period as extracted zg file): ', np.shape(pr))
print ('##################XXXXX##################')


# Indices that correspond to 29 Feb of leap years within 1980-2005 period
# gdel for year Y computed as the sum of:
#      1. 365 X no.s of non-leap years before Y
#      2. 366 X no.s of leap years before Y
#      3. 31+28 (Jan and Feb months of Y)

# Removing 29 Feb from Leap years
if g == 1:
    gdel = [59, 1520, 2981, 4442, 5903, 7364, 8825] # indices for 29 feb of leap years
    pr = np.delete(pr, gdel, 0)
    print ('size of extracted prec file after removing 29 febs: ', np.shape(pr))        

   
################################################################################
################################################################################
#### zg time slice ####
print ('##### extract GEOPOTENTIAL ######')
# number of lag days
pl = 4
zst = '250'

# extract lat range and dim
flnamez = cpathz+"/"+flz[0]
ncz = NetCDFFile(flnamez, 'r')

plev = ncz.variables['plev'][pl]
print ('EXTRACTING DATA FOR LEVEL :: ', plev, ' Pa')

lat1 = ncz.variables['lat'][:]
lonz = ncz.variables['lon'][:]
ilt = np.where ((lat1>=-90) & (lat1<=20))
latz = lat1[ilt[0]]
zgtest = ncz.variables['zg'][0,pl,ilt[0],:]
[nz,oz] = np.shape(zgtest)
del zgtest

zg = np.zeros((1,nz,oz))
gdel = []; g1 = 0; g2 = 0;

for ip in range(len(flz)): # loop over all the zg files
    flpr = cpathz+"/"+flz[ip]
    print (flz[ip])
    nc = NetCDFFile(flpr, 'r')
    tim = nc.variables['time']
    tim = tim[:] - time_offset
    tim1 = tim[0]
    tim2 = tim[-1]
    print ('1st & last time steps:: ', tim1, tim2)
    
# find the time slices we need from each file

    if tim2 < t1:
        continue # go to next iteration
    elif tim1<=t1:
        int1 = np.where(tim[:] == t1)
        int1 = int1[0][0] # tuple to integer
    elif tim1>t1:
        int1 = 0

    if tim2 <= t2-1:  # t2-1 == 31 Dec of 30th year
        int2 = len(tim)-1
    elif tim2 > t2-1:
        int2 = np.where(tim[:] == t2)
        int2 = int2[0][0] # tuple to integer
    
    print ('1st & last indices from this file :: ', int1,int2)
    pr1 = nc.variables['zg'][int1:int2+1,pl,ilt[0],:] # timeslice up to (t2-1)//last index is intm in actuality i.e., a[0:5] will give a[0], a[1], a[2], a[3], a[4] i.e., 5 elements
    print ('size of extracted zg from this file ::', np.shape(pr1))
    
       
    zg = np.vstack((zg,pr1))
    del pr1
    
    if tim2 >= t2-1:
        break

zg = zg[1:,:,:]
print ('size of extracted zg file (for same period as extracted prec file): ', np.shape(zg))
print ('##################XXXXX##################')


# Indices that correspond to 29 Feb of leap years within 1980-2005 period
# gdel for year Y computed as the sum of:
#      1. 365 X no.s of non-leap years before Y
#      2. 366 X no.s of leap years before Y
#      3. 31+28 (Jan and Feb months of Y)

# gdel
if g == 1:
    gdel = [59, 1520, 2981, 4442, 5903, 7364, 8825] # indices for 29 feb of leap years
    zg = np.delete(zg, gdel, 0)
    print ('size of extracted zg file after removing 29 febs: ', np.shape(zg))        


##################### SEGMENT - 3 ####################
##### finding anomaly/climatology from the extracted datasets
######################################################

[m,n,o] = np.shape(zg)

#time = np.linspace(1,3600,3600)
time = np.linspace(1,m,m)

print ('#### Computing daily climatology ####')
# computing daily climatology
zg_clim1 = np.zeros((int(dy),nz,oz))
pr_clim1 = np.zeros((int(dy),npr,opr))


tsp = int(dy)
for iyr in range(int(len(zg)/dy)):
    yrsp = iyr*int(dy)
    ist, iend = yrsp, yrsp + tsp
    zgm = zg[ist:iend,:,:]
    prm = pr[ist:iend,:,:]
    zg_clim1 = np.nansum((zg_clim1,zgm), axis=0)
    pr_clim1 = np.nansum((pr_clim1,prm), axis=0)
    time1 = time [ist:iend]
    print ('start and end of each model year : ', time1[0], time1[-1])
    print ('length of each year : ', time1[-1] - time1[0]+1)

zg_clim = zg_clim1/(iyr+1)
pr_clim = pr_clim1/(iyr+1)

del zg_clim1
del pr_clim1
print ('shape of daily zg annual cycle:: ', np.shape(zg_clim))
print ('max value of zg clim: ', np.nanmax(zg_clim[:]))
print ('shape of daily pr annual cycle:: ', np.shape(pr_clim))
print ('max value of pr clim: ', np.nanmax(pr_clim[:]))

# * * * * * * * * * * * * * * * * *
# WITH FFT SMOOTHING!!
# using rfft: to smooth climatology
# https://stackoverflow.com/questions/23077850/fourier-smoothing-of-data-set
# https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.rfft.html#numpy.fft.rfft

############ zg
zfft = np.fft.rfft(zg_clim, axis=0)
zfft[6:,:,:] = 0
zg_smooth1 = np.fft.irfft(zfft, axis=0, n=int(dy)) # specify the length of final output to 365/360 *** very imp for odd no. samples (i.e., 365)
print (np.shape(zg_smooth1))
############ pr
pfft = np.fft.rfft(pr_clim, axis=0)
pfft[6:,:,:] = 0
pr_smooth1 = np.fft.irfft(pfft, axis=0, n=int(dy)) # specify the length of final output to 365/360 *** very imp for odd no. samples (i.e., 365)
print (np.shape(pr_smooth1))

print ('shape of zg smooth: ', np.shape(zg_smooth1))
print ('shape of pr smooth: ', np.shape(pr_smooth1))

#del zg_clim; del pr_clim


#plt.close("all")
#fig = plt.figure
#ax1=plt.subplot(211)
#ax1.plot(zg_smooth1[:,23,32], 'k')
#ax1.plot(zg_clim[:,23,32], 'b')
#plt.legend(fontsize=10)
#figname = ('test_fft.png')
#plt.savefig(figname, dpi = 600)

# repeat the annual cycle 10 times (= no.s of years) so that
# can be readily subtracted from prc (daiy data) to compute daily anomaly

################### SEGMENT - 4 ########################
### finding anomaly
# ZG
# tiling clim data
zg_smooth = np.tile(zg_smooth1,(int(len(zg)/dy),1,1))
del zg_smooth1
# finding anomaly
zgano1 = zg - zg_smooth
del zg; del zg_smooth
print ('shape of daily anomaly var - for all 365/360 days x 20 yearsi:: ', np.shape(zgano1))




# import sys 
# sys.exit()




# removing linear trend from anomaly data
#zgano=scipy.signal.detrend(zgano1, axis=0)
#print 'shape of detrended anomaly var - for all 360 days x 30 years'
#print np.shape(zgano)

# save zgano, then delete :: for memory issues
# saving the unfiltered anomaly data (comment if not interested in raw anomaly data)
flano = ('zg'+zst+'_ano_daily_detrend_NoFilt_SH_'+modl+'_'+fyr+'.npz')
print (' raw/unfiltered anomaly data saved in::: ', flano)
np.savez(flano, zgano=zgano1, time=time, lon=lonz, lat=latz)
del zgano1


# PREC
# tiling clim data
pr_smooth = np.tile(pr_smooth1,(int(len(pr)/dy),1,1))
print ('shape of annual cycle var, repeated 10 times:: ', np.shape(pr_smooth))
del pr_smooth1
# finding anomaly
prano1 = pr - pr_smooth
del pr; del pr_smooth
print ('shape of daily anomaly var - for all 365/360 days x 20 yearsi:: ', np.shape(prano1))

#########################################################

# daily anomaly for all timesteps (i.e., for all 10 years)

# removing linear trend from anomaly data
#prano=scipy.signal.detrend(prano1, axis=0)
#print 'shape of detrended anomaly var - for all 365/360 days x 20 years'
#print np.shape(prano)



# save prano, then delete :: for memory issues
# saving the unfiltered anomaly data (comment if not interested in raw anomaly data)
flano = ('pr_ano_daily_detrend_NoFilt_NH_'+modl+'_'+fyr+'.npz')
print (' raw/unfiltered anomaly data saved in::: ', flano)
np.savez(flano, prano=prano1, time=time, lon=lonp, lat=latp)
# del prano1



