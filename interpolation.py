# importing required libraries
import pandas as pd
import matplotlib as mlp
from matplotlib import pyplot as plt
import numpy as np
#import seaborn as sbd
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import xarray as xr
from shapely.geometry import Point
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import geopandas as gpd
import geodatasets

plt.close()

fileLoc = '/g100_work/OGS23_PRACE_IT/gbuccino/test/Venice_Lagoon_hydrodynamics/'
area = 'DATA_Venice_lagoon' #can be Venice_Lagoon, Grado_Lagoon 

bathy_grid = pd.read_csv(f'{fileLoc}'+f'{area}.csv', skiprows=0, header=0)
bathy_unstruct = pd.read_csv(f'{fileLoc}'+'bathy_2003CORILA_4326.csv',names= ['idx','lon','lat','depth1','depth2', 'geometry'],header=0)
#bathy_chan = pd.read_csv(f'{fileLoc}'+'bathy_coarsed2013_4326.csv',skiprows=0)

# Initialization
latbins = bathy_grid['lat'][0:111:1]
lonbins = bathy_grid['lon'][0::111]


lons_un = np.load(f'{fileLoc}'+'lon_2003CORILA_4326.npy')        #unstructured data
lats_un = np.load(f'{fileLoc}'+'lat_2003CORILA_4326.npy')        #unstructured data
depth_un = np.load(f'{fileLoc}'+'depth_2003CORILA_4326.npy') 
# lons_chan = bathy_chan['geometry'].x        #channel data
# lats_chan = bathy_chan['geometry'].y        #channel data

sum_temp = np.zeros((latbins.size-1, lonbins.size-1))
occ_temp = np.zeros((latbins.size-1, lonbins.size-1))

occ_temp, lonbins, latbins = np.histogram2d(lons_un, lats_un, bins = [lonbins,latbins])
sum_temp, lonbins, latbins = np.histogram2d(lons_un, lats_un, bins = [lonbins, latbins], weights=depth_un)
sum_temp[sum_temp>0] = None


occ_temp = occ_temp.T
sum_temp = sum_temp.T

# Compute mean
#FILL_VALUE = -9999
mean_temp = np.ones((latbins.size-1, lonbins.size-1))*2
mean_temp[occ_temp != 0] = sum_temp[occ_temp != 0]/occ_temp[occ_temp != 0]

np.save(f'{fileLoc}'+f'bathy_interpolated_dataBin.npy', mean_temp)

mean_temp = np.ma.masked_where(mean_temp >=0, mean_temp)
occ_temp= np.ma.masked_where(occ_temp>=0, occ_temp)
#mean_temp[mean_temp == 0] = None

# ****************
# *** Plot map ***
X, Y = np.meshgrid(lonbins, latbins)
plt.figure(figsize=(5, 10))

# Plot mean data
ax = plt.subplot(211)
ax.set_facecolor('0.8')
plt.pcolormesh(X, Y, mean_temp, cmap=cm.RdBu_r)
plt.colorbar(label='mean depth into the pixel')
plt.clim(np.min(sum_temp), 0.0000)
# plt.xlim(1, 4)
# plt.ylim(40, 44)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Plot number of observations

ax = plt.subplot(212)
plt.pcolormesh(X, Y, occ_temp)
plt.colorbar(label='Number of observations')
# plt.xlim(1, 4)
# plt.ylim(40, 44)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.savefig(f'{fileLoc}'+'mean_bathy.pdf', format='pdf')
plt.close()





print('ciao')