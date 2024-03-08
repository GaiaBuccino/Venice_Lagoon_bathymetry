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

def geom_creation(coord1: np.ndarray, coord2:np.ndarray):
    points = []

    for ii in coord1:
        for jj in coord2:
            points.append(Point(ii,jj))

    return points

x_dim = 111
y_dim = 111

fileLoc = '/g100_work/OGS23_PRACE_IT/gbuccino/test/'
area = 'Venice_lagoon_DATA' #can be Venice_Lagoon, Grado_Lagoon 
fileName = fileLoc + f'bathy_{area}.bin'
fileLat = fileLoc + f'lat_{area}.bin'
fileLon = fileLoc + f'lon_{area}.bin'

data_bat = np.fromfile(fileName, dtype="<f4", count=-1).reshape(x_dim, y_dim)
data_lon = np.fromfile(fileLon, dtype="<f4", count=-1)
data_lat = np.fromfile(fileLat, dtype="<f4", count=-1)

da_bath = xr.DataArray(name='depth', data =data_bat, dims = ("lat","lon"), coords={"lat": data_lat, "lon": data_lon},)


data_pd = da_bath.to_dataframe()
data_csv = data_pd.to_csv(f'{fileLoc}'+f'{area}_Querin.csv')

# coords_Q = dict({'lat':data_lat, 'lon':data_lon})
# map_Q = np.meshgrid(data_lon, data_lat) 

""" plt.pcolormesh(data_bat)
plt.colorbar()
plt.clim(0,-4)
plt.show()

plt.savefig(f'{fileLoc}Bathy_{area}') """


## INTERPOLATION

#contructing geometry:
pairs = (data_lat,data_lon)
data_pd.insert(0, 'geometry', pairs)
bathy_Q = gpd.GeoDataFrame(data_pd, geometry = 'geometry', crs= "EPSG:4326")


bathy_C_pd = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/Bati_2013_coarsed.csv',header=None,sep=',')
bathy_C_pd.insert(0, 'geometry', geom_creation(bathy_C_pd[0], bathy_C_pd[1]))
bathy_C = gpd.GeoDataFrame(bathy_C_pd, geometry = 'geometry', crs= "EPSG:3004")
bathy_C_4326 = bathy_C.to_crs("EPSG:4326")
bathy_C_4326.to_csv('bathy_2013_4326.csv')
#bathy_C_chan = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/Batimetria_2003CORILA_shift_40-15.csv',header=None,sep=',')


