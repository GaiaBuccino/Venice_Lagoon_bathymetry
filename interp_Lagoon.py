#import struct
import pandas as pd
import matplotlib as mlp
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sbd
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import xarray as xr
import geopandas

 
# importing required libraries
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


bathy_Q = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/ONLY_Venice_lagoon_Querin.csv',header=None,sep=',')
bathy_C = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/Bati_2013_coarsed.csv',header=None,sep=',')
bathy_C_chan = pd.read_csv('/g100_work/OGS23_PRACE_IT/gbuccino/test/Batimetria_2003CORILA_shift_40-15.csv',header=None,sep=',')

da_bath_Q = xr.DataArray(name='depth', data =bathy, dims = ("lat","lon"), coords={"lat": data_lat, "lon": data_lon},)


DATAFRAME_crs3004= geopandas.GeoDataFrame(bathy_Q, geometry='coords',crs=3004)
DATAFRAME_crs3035 =DATAFRAME_crs3004.to_crs(3004)