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
bathy_new = bathy_grid[['lat', 'lon']].copy()


# Initialization
lat_file = np.array(bathy_grid['lat'][0::111])
lon_file = np.array(bathy_grid['lon'][0:111:1])

# Enlarge the grid
dlon = lon_file[1] - lon_file[0]
dlat = lat_file[1] - lat_file[0]

lonb = lon_file - 0.5*dlon
latb = lat_file - 0.5*dlat

lonbins = np.zeros((len(lonb)+1))
lonbins[0:len(lonb)] = lonb
lonbins[-1] = np.max(lon_file) + 0.5*dlon

latbins = np.zeros((len(latb)+1))
latbins[0:len(latb)] = latb
latbins[-1] = np.max(lat_file) + 0.5*dlat


xx, yy = np.meshgrid(latbins,lonbins)

# Upload files

lons_un = np.load(f'{fileLoc}'+'lon_2003CORILA_4326.npy')        #unstructured data
lats_un = np.load(f'{fileLoc}'+'lat_2003CORILA_4326.npy')        #unstructured data
depth_un = np.load(f'{fileLoc}'+'depth_2003CORILA_4326.npy') 
# lons_chan = bathy_chan['geometry'].x        #channel data
# lats_chan = bathy_chan['geometry'].y        #channel data

sum_dep = np.zeros((latbins.size-1, lonbins.size-1))
occ_dep = np.zeros((latbins.size-1, lonbins.size-1))

occ_dep, lonb, latbins = np.histogram2d(lons_un, lats_un, bins = [lonbins,latbins])
sum_dep, lonb, latbins = np.histogram2d(lons_un, lats_un, bins = [lonbins, latbins], weights=depth_un)
#sum_dep[sum_dep>0] = None    fare dopo il conto della media

occ_dep = occ_dep.T
sum_dep = sum_dep.T

# Compute mean
#FILL_VALUE = -9999
mean_dep = np.ones((latbins.size-1, lonbins.size-1))*10
mean_dep[occ_dep != 0] = sum_dep[occ_dep != 0]/occ_dep[occ_dep != 0]

min_avg_dep = np.min(mean_dep)

mean_dep[occ_dep == 0] = None
mean_dep[mean_dep>0] = None

dep = mean_dep.flatten()

bathy_new['depth'] = pd.Series(dep) 
bathy_grid[bathy_new['depth'].isna()] = None


""" # Choose a pixel (minimum in avg depth)

check =  min_avg_dep
x_check, y_check = np.where(mean_dep == check)

hist_check_x = [lonbins[x_check[0]], lonbins[x_check[0]+1]]
hist_check_y = [latbins[y_check[0]], latbins[y_check[0]+1]]

occ_check, lo, la = np.histogram2d(lons_un, lats_un, bins = 1, range = [hist_check_x, hist_check_y])
sum_check, lo, la = np.histogram2d(lons_un, lats_un, bins = 1, range = [hist_check_x, hist_check_y], weights=depth_un)

binEdges, np.digitize(A, binEdges)


x,y = np.meshgrid(hist_check_x, hist_check_y)
# plt.xlim(lonbins[x_check[0]], lonbins[x_check[0]+1])
# plt.ylim(latbins[y_check[0]], latbins[y_check[0]+1])
plt.grid()
plt.plot(x,y, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")
#plt.plot(occ_check[], marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")
plt.savefig(f'{fileLoc}'+'prova_dots.pdf', format='pdf')
plt.close()
 """


""" new_bathy_df = { 'lons': xx,
                 'lats': yy,
                 'depth': mean_dep

}

new_bathy_df = pd.DataFrame(new_bathy_df)
new_bathy_df.to_csv(f'{fileLoc}'+f'bathy_interpolated_dataBin.csv')
 """
# mean_dep = np.ma.masked_where(occ_dep==0, mean_dep)
# occ_dep= np.ma.masked_where(occ_dep==0, occ_dep)

#bathy_new['depth'] = 


# ****************
# *** Plot map ***
X, Y = np.meshgrid(lonbins, latbins)
plt.figure(figsize=(5, 10))

# Plot mean data
ax = plt.subplot(211)
ax.set_facecolor('0.8')
plt.pcolormesh(X, Y, mean_dep, cmap=cm.RdBu_r)
plt.colorbar(label='mean depth into the pixel')
plt.clim(np.min(mean_dep), 0.0000)
# plt.xlim(1, 4)
# plt.ylim(40, 44)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Plot number of observations

ax = plt.subplot(212)

plt.pcolormesh(X, Y, occ_dep)
plt.colorbar(label='Number of observations')
# plt.xlim(1, 4)
# plt.ylim(40, 44)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.savefig(f'{fileLoc}'+'mean_bathy.pdf', format='pdf')
plt.close()

# Plot data on subplots
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax1.title.set_text('Interpolated bathymetry on DATA')
plt.pcolormesh(X,Y,bathy_new['depth'].to_numpy().reshape(111,111))
plt.colorbar()
plt.clim(min_avg_dep, 0.0)
plt.savefig(f'{fileLoc}'+ 'Interpolated_on_DATAbin.pdf')


ax2 = fig.add_subplot(132)
ax2.title.set_text('Original bathymetry')
plt.pcolormesh(X,Y,bathy_grid['depth'].to_numpy().reshape(111,111))
plt.colorbar()
plt.clim(min_avg_dep, 0.0)
plt.savefig(f'{fileLoc}'+ 'Original_DATAbin.pdf')


bathy_diff = bathy_grid[['lat', 'lon','depth']].copy()
bathy_diff['depth'] = (bathy_grid['depth'] - bathy_new['depth'])
ax3 = fig.add_subplot(133)
ax3.title.set_text('Difference in bathymetry')
plt.pcolormesh(X,Y,bathy_diff['depth'].to_numpy().reshape(111,111))
plt.colorbar()
plt.savefig(f'{fileLoc}'+ 'Bathymetries_difference.pdf', format='pdf')

fig.savefig(f'{fileLoc}'+ 'Bathymetries_comparison.pdf', format='pdf')


""" plt.figure(figsize=(30, 30))
ax1 = fig.add_subplot(131)
plt.pcolormesh(X,Y,bathy_new['depth'].to_numpy().reshape(111,111))
plt.colorbar()
plt.clim(min_avg_dep, 0.0)
ax1.title.set_text('Interpolated bathymetry on DATA')
plt.savefig(f'{fileLoc}'+ 'Bathymetry_Interpolated_on_DATA.pdf', format='pdf')

ax2 = fig.add_subplot(132)
plt.pcolormesh(bathy_grid['lon'].to_numpy().reshape(111,111),bathy_grid['lat'].to_numpy().reshape(111,111),bathy_grid['depth'].to_numpy().reshape(111,111))
plt.colorbar()
ax2.title.set_text('Original bathymetry')

ax3 = fig.add_subplot(133)
bathy_diff = bathy_grid[['lat', 'lon','depth']].copy()
bathy_diff['depth'] = (bathy_grid['depth'] - bathy_new['depth'])
plt.pcolormesh(bathy_new['lon'].to_numpy().reshape(111,111),bathy_new['lat'].to_numpy().reshape(111,111),bathy_diff['depth'].to_numpy().reshape(111,111))
plt.colorbar()
ax3.title.set_text('Difference among bathymetry')

plt.savefig(f'{fileLoc}'+ 'Bathymetry_on_DATA.pdf', format='pdf') """


plt.close()



print('ciao')