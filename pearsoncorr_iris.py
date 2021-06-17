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
import iris
from zg2psi import zg2psi
from iris.analysis.stats import pearsonr
from iris.analysis.cartography import area_weights
#######################################################
    
cubes = iris.load('test_files_1.nc')
zg = cubes[0]
print(zg)
 
wg = iris.analysis.cartography.area_weights(zg, normalize=False)
 
kr = iris.analysis.stats.pearsonr (zg,zg, corr_coords=None, weights=None, mdtol=1.0, common_mask=False)