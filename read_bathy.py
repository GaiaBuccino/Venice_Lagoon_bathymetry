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
import time

def geom_creation(coord1: np.ndarray, coord2:np.ndarray, structured:bool = False): #-> np.ndarray[Point]
    
    if structured:
        points = []
        for jj in coord2:
            for ii in coord1:

                pt = Point(jj,ii)    #coordinates of points given as lat, lon
                points.append(pt)

        return points

    else:
        pandaspoints = pd.DataFrame({'lon':coord2, 'lat':coord1})
        pandaspoints['coords'] = list(zip(pandaspoints['lon'],pandaspoints['lat']))
        pandaspoints['coords'] = pandaspoints['coords'].apply(Point)
        return pandaspoints['coords']


####  PREPARATION OF DATA & SAVING

## STRUCTURED GRID DATA
    
x_dim = 792   #from DATA.bin 111 x 111, from RESCALED.bin 61x61, 792x424 from RESCALED1_128_NO_lagoon, 61x61 from _cut_Venice_lagoon
y_dim = 424

fileLoc = '/g100_work/OGS23_PRACE_IT/gbuccino/test/Venice_Lagoon_hydrodynamics/'
area = 'completa_1_128' #can be RESCALED1_128_  OR  DATA_ OR lagoon_ADRIsize OR ADRI_1_128_NO_Lagoon


fileName = fileLoc + f'bathy_{area}.bin' 
fileLat = fileLoc + f'lat_{area}.bin'
fileLon = fileLoc + f'lon_{area}.bin'   

data_bat = np.fromfile(fileName, dtype="<f4", count=-1).reshape(y_dim, x_dim)
data_lon = np.fromfile(fileLon, dtype="<f4", count=-1)
data_lat = np.fromfile(fileLat, dtype="<f4", count=-1)

da_bath = xr.DataArray(name='depth', data =data_bat, dims = ("lat","lon"), coords={"lon": data_lon,"lat": data_lat})  # coords have to be in the same order of the file data_bat

data_pd = da_bath.to_dataframe()
data_pd.to_csv(f'{fileLoc}'+f'{area}.csv')
# bathy_Q = gpd.GeoDataFrame(data_pd, geometry=geom_creation(data_lon,data_lat,structured = True), crs="EPSG:4326") 
# bathy_Q.to_csv(f'{fileLoc}'+f'{area}_geodata.csv')

plt.pcolormesh(data_bat)
plt.colorbar()
plt.clim(0,-4)
plt.show()

plt.savefig(f'{fileLoc}Bathy_{area}') 

## UN-STRUCTURED DATA

bathy_C_pd = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/Venice_Lagoon_hydrodynamics/Batimetria_2003CORILA_shift_40-15.csv',header=None,sep=',')
bathy_C = gpd.GeoDataFrame(bathy_C_pd, geometry = geom_creation(bathy_C_pd[1], bathy_C_pd[0]), crs= "EPSG:3004")
# bathy_C_lat3004 = bathy_C['geometry'].y
# bathy_C_lon3004 = bathy_C['geometry'].x
#bathy_C.to_csv(f'{fileLoc}'+f'bathy_2003CORILA_3004.csv')

bathy_C_4326 = bathy_C.to_crs("EPSG:4326")
bathy_C_lat4326 = np.array(bathy_C_4326['geometry'].y)
bathy_C_lon4326 = np.array(bathy_C_4326['geometry'].x)
bathy_C_depth4326 = np.array(bathy_C_4326[2])           # depth presa dalla seconda colonna 
bathy_C_4326.to_csv(f'{fileLoc}'+f'bathy_2003CORILA_4326.csv')

np.save(f'{fileLoc}'+f'lat_2003CORILA_4326.npy', bathy_C_lat4326)
np.save(f'{fileLoc}'+f'lon_2003CORILA_4326.npy',bathy_C_lon4326)
np.save(f'{fileLoc}'+f'depth_2003CORILA_4326.npy',bathy_C_depth4326)

bathy_C_chan_pd = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/Venice_Lagoon_hydrodynamics/Bati_2013_coarsed.csv',header=None,sep=',')
bathy_C_chan = gpd.GeoDataFrame(bathy_C_chan_pd, geometry = geom_creation(bathy_C_chan_pd[1], bathy_C_chan_pd[0]), crs= "EPSG:3004")
bathy_C_chan.to_csv(f'{fileLoc}'+f'bathy_coarsed2013_3004.csv')
bathy_C_chan_4326 = bathy_C_chan.to_crs("EPSG:4326")

bathy_C_chan_lat4326 = bathy_C_chan_4326['geometry'].y
bathy_C_chan_lon4326 = bathy_C_chan_4326['geometry'].x
bathy_C_chan_depth4326 = np.array(bathy_C_chan_4326[2])

np.save(f'{fileLoc}'+f'lat_coarsed2013_4326.npy',bathy_C_lat4326)
np.save(f'{fileLoc}'+f'lon_coarsed2013_4326.npy',bathy_C_lon4326)
np.save(f'{fileLoc}'+f'depth_coarsed2013_4326.npy',bathy_C_chan_depth4326)

bathy_C_chan_4326.to_csv(f'{fileLoc}'+f'bathy_coarsed2013_4326.csv')
