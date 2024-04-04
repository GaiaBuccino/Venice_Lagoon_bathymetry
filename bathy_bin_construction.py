# importing required libraries
import pandas as pd
import matplotlib as mlp
from matplotlib import pyplot as plt
import numpy as np
#import seaborn as sbd
from matplotlib import cm
# from matplotlib.ticker import LinearLocator
# import xarray as xr
# from shapely.geometry import Point
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


plt.close()

fileLoc = '/g100_work/OGS23_PRACE_IT/gbuccino/test/Venice_Lagoon_hydrodynamics/'
#area = 'RESCALED1_128_Venice_lagoon' 

x_dim = 792   #from DATA.bin 111 x 111, from RESCALED.bin 61x61, 792x424 from RESCALED1_128_NO_lagoon, 61x61 from _cut_Venice_lagoon
y_dim = 424

bathy_Adri = pd.read_csv(f'{fileLoc}'+'ADRI_1_128_NO_Lagoon.csv', skiprows=0, header=0)
lat_Adri = np.array(bathy_Adri['lat'][0:335807:792])
bathy_Adri[bathy_Adri['depth'].isna()] = 0.0
lon_Adri = np.array(bathy_Adri['lon'][0:792])
#dep_Adri = np.array(bathy_Adri['depth']).reshape(x_dim,y_dim)


bathy_temp = pd.read_csv(f'{fileLoc}'+'adjusted_bathy.csv', skiprows=0, header=0)
bathy_temp[bathy_temp['depth'].isna()] = 0.0
bathy_temp.to_csv(f'{fileLoc}' + 'Final_bathymetry_only_Lagoon.csv', index=False)

#bathy_temp[bathy_adj['depth'].isna()] = 0.0

bathy_adj = bathy_Adri[['lat', 'lon', 'depth']].copy()
bathy_adj['depth'] = bathy_temp['depth'] + bathy_Adri['depth']
bathy_adj.to_csv(f'{fileLoc}' + 'Final_bathymetry.csv', index=False)



bathy_completa = pd.read_csv(f'{fileLoc}'+'completa_1_128.csv', skiprows=0, header=0)
# lon_adj = lon_Adri
# lat_adj = lat_Adri
# dep_adj = np.array(bathy_adj['depth'])

# bathy_final = bathy_Adri[['lat', 'lon']].copy()
# bathy_final['depth'] = bathy_adj['depth'] + bathy_Adri['depth']


xx, yy = np.meshgrid(lon_Adri,lat_Adri)

fig = plt.figure()
plt.pcolormesh(xx,yy,bathy_temp['depth'].to_numpy().reshape(y_dim,x_dim))
plt.colorbar()
plt.clim(-4, 0.0000)
plt.savefig(f'{fileLoc}'+'Completa_only_lagoon_bathymetry')
plt.close()


# fig2 = plt.figure()
# plt.pcolormesh(xx,yy,bathy_completa['depth'].to_numpy().reshape(y_dim,x_dim))
# plt.colorbar()
# plt.clim(-4, 0.0000)
# plt.savefig(f'{fileLoc}'+'Original_1_128_bathymetry')
# plt.close()

# fig3 = plt.figure()
# plt.pcolormesh(xx,yy,(bathy_completa['depth'].to_numpy().reshape(y_dim,x_dim)-bathy_adj['depth'].to_numpy().reshape(y_dim,x_dim)))
# plt.colorbar()
# plt.clim(-4, 0.0000)
# plt.savefig(f'{fileLoc}'+'Difference_completa-adj_1_128_bathymetry')
# plt.close()

print('ciau')


