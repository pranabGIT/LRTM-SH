import iris
import numpy as np

def zg2psi(zgCube):
    "Convert from geopotential height into streamfunction"
    # Create an iris cube of the coriolis frequency
    latCoord = zgCube.coord('latitude').points
    lonCoord = zgCube.coord('longitude').points
    lon2D,lat2D  = np.meshgrid(lonCoord,latCoord)
    omega = 7.29e-5 #Rotation frequency
    # Corolis paramter
    f0 = 2.*omega*np.sin(np.deg2rad(lat2D))
    #create_cube
    coriolisCube = create_cube(f0,zgCube)
    coriolisCube.standard_name = 'coriolis_parameter'
    coriolisCube.units = 's-1'
    # Convert \psi = zg/f_0
    psiCube = iris.analysis.maths.divide(zgCube,coriolisCube)
    return psiCube
