# importing required libraries
from math import nan
import pandas as pd
#import matplotlib as mlp
from matplotlib import pyplot as plt
import numpy as np
# import seaborn as sbd
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator
import xarray as xr
from shapely.geometry import Point
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import geopandas as gpd
import os
#import time
from typing import List


# def csv_creation(fileLoc: str, area: str, fileDest: str, x_dim: int, y_dim: int, plot:bool = False):
    
#     """Function that creates the csv files to store data ([lon, lat depth]) starting from files bin

#     Args:
#         fileLoc (str): path of the file
#         area (str): area under analysis (required 'bathy_{area}.bin', 'lat_{area}.bin', 'lon_{area}.bin')  
#         fileDest (str): path where the produced fles have to be saved
#         x_dim (int): points of the grid in the x direction
#         y_dim (int): points of the grid in the y direction
#         plot (bool, optional): variable representing if a plot has to be produced and saved. Defaults to False.
#     """
    

#     fileName = fileLoc + f'bathy_{area}.bin' 
#     #fileLat = fileLoc + f'lat_{area}.bin'
#     # fileLon = fileLoc + f'lon_{area}.bin'   

    

#     data_bat = np.fromfile(fileName, dtype="<f4", count=-1).reshape(x_dim, y_dim)
#     data_lon = np.fromfile(fileLon, dtype="<f4", count=-1)
#     data_lat = np.fromfile(fileLat, dtype="<f4", count=-1)

#     da_bath = xr.DataArray(name='depth', data =data_bat, dims = ("lat","lon"), coords={"lon": data_lon,"lat": data_lat})  # coords have to be in the same order of the file data_bat

#     data_pd = da_bath.to_dataframe()
#     data_pd.to_csv( f'{fileDest}'+f'{area}.csv')
    
#     if plot:
#         plt.pcolormesh(data_bat)
#         plt.colorbar()
#         plt.clim(0,-4)
#         plt.show()

#         plt.savefig(f'{fileDest}Bathy_{area}') 

#     return

    
def geom_creation(coord1: np.ndarray, coord2:np.ndarray, structured:bool = False): #-> np.ndarray[Point]

    """Function that creates the structure "geometry" required to change the reference system of coordinates

    Args:
        coord1 (np.ndarray): longitudes
        coord2 (np.ndarray): latitudes
        structured (bool, optional): variable representing if the grid is structured or not. Defaults to False -> unstructured.

    Returns:
        points: array of Point representing the geometry
    """

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



def change_ref_system(bathy_pd: pd.DataFrame, fileDest: str , fileName:str, initial_ref: str, final_ref: str):
    """Function that converts a pandas DataFrame into a GeoDataFrame structure (it adds the geometry variable and allows changing of the reference system in an easy way)

    Args:
        file_path (str): path where the file is stored
        fileDest (str): path where the created files have to be saved
        fileName (str): name for the produced filesgbu
        initial_ref (str): starting system of coordinates 
        final_ref (str): final system of coordinates (suppose a 'EPSG:' before)
    """
    converted_path = fileDest +f'{fileName}_{final_ref}.csv'
    
    if not os.path.exists(converted_path):
        
        coord_start = "EPSG:"+ initial_ref
        geom = geom_creation(bathy_pd[1], bathy_pd[0])
        bathy_gpd = gpd.GeoDataFrame(bathy_pd, geometry = geom, crs = coord_start ) # type: ignore

        bathy_gpd_converted = bathy_gpd.to_crs(final_ref)
        lat_converted = np.array(bathy_gpd_converted['geometry'].y)
        lon_converted = np.array(bathy_gpd_converted['geometry'].x)
        
        if fileName == 'Bathy_2003CORILA':
            bathy_gpd_converted = bathy_gpd_converted.drop(bathy_gpd_converted.columns[2], axis =1)
            bathy_gpd_converted = bathy_gpd_converted.drop(bathy_gpd_converted.columns[1], axis =1)
            bathy_gpd_converted = bathy_gpd_converted.drop(bathy_gpd_converted.columns[0], axis =1)
            bathy_gpd_converted.insert(0, "lon", lon_converted)
            bathy_gpd_converted.insert(1, "lat", lat_converted)
            bathy_gpd_converted.rename(columns={0: "lon", 1: "lat", 3:"depth", "geometry": "geometry"}, inplace = True)
            #bathy_gpd_final['depth'].loc[bathy_gpd_final.depth > 0] = 0.0
            bathy_gpd_converted.loc[bathy_gpd_converted["depth"] > 0.0, "depth"] = 0.0
            #bathy_gpd_final[bathy_gpd_final.depth > 0] = 0.0
            bathy_gpd_converted.to_csv(fileDest +f'{fileName}_{final_ref}.csv', header= ['lon','lat','depth', 'geometry'], index=False)

        elif fileName == 'Bathy_2013_coarsed':
            bathy_gpd_converted = bathy_gpd_converted.drop(bathy_gpd_converted.columns[1], axis =1)
            bathy_gpd_converted = bathy_gpd_converted.drop(bathy_gpd_converted.columns[0], axis =1)
            bathy_gpd_converted.insert(0, "lon", lon_converted)
            bathy_gpd_converted.insert(1, "lat", lat_converted)
            bathy_gpd_converted.rename(columns={0: "lon", 1: "lat", 2:"depth", "geometry": "geometry"}, inplace = True)
            bathy_gpd_converted.loc[bathy_gpd_converted["depth"] > 0.0, "depth"] = 0.0    #per sicurezza
            bathy_gpd_converted.to_csv(fileDest +f'{fileName}_{final_ref}.csv', header= ['lon','lat','depth', 'geometry'], index=False)

    return f'{fileName}_{final_ref}'        #lat_fref, lon_fref, depths, bathy_gpd_final


def interp_on_structured(lon_st:np.ndarray, lat_st:np.ndarray, files: List[pd.DataFrame], fileName: List[str], nPt_int: int, fileDest: str, perc: int):
    """Given a list of files with unstructured data (csv), thay are interpolated on a structured grid through the computation of the spatial average on a more refined 
    grid in order to have more accurate values

    Args:
        lon_st (np.ndarray): array with the longitude data of the structured grid on which the interpolation has to be performed
        lat_st (np.ndarray): array with the latitude data of the structured grid on which the interpolation has to be performed
        files (List[str]): list of the name of the file with the unstructured data tht have to be interpolated on the structured grid
        nPt_int (int): number of points on each cell to perform the refinement of the grid (RM: nPt points means nPt-1 refinements)
        fileDest (str): string with the destination folder where the produced files are saved

    Returns:
        xr.DataArray: file with the grid values (lon, lat) and the averaged depths
    """

    nRef = nPt_int-1

    depth_path = fileDest + f'Averaged_depth_{nRef}ref.nc'
    percentage_path = fileDest + f'Water_coverage_{perc}perc_{nRef}ref.csv'

    if not (os.path.exists(depth_path) & os.path.exists(percentage_path)):

        # Enlarge the grid to have values in the center of the cells instead of on their vertices
        dlon = lon_st[1] - lon_st[0]
        dlat = lat_st[1] - lat_st[0]

        lonb = lon_st - 0.5*dlon
        latb = lat_st - 0.5*dlat

        lonbins = np.zeros((len(lonb)+1))
        lonbins[0:len(lonb)] = lonb
        lonbins[-1] = np.max(lon_st) + 0.5*dlon
        latbins = np.zeros((len(latb)+1))
        latbins[0:len(latb)] = latb
        latbins[-1] = np.max(lat_st) + 0.5*dlat

        # Definition of the matrices required for the computation

        # Global matrices, same dimension of the structured grid (len(latbins)-1,len(lonbins)-1)
        # tot_occ = np.zeros((len(latbins)-1,len(lonbins)-1))
        # perc_occ = np.zeros((len(latbins)-1,len(lonbins)-1))
        global_percentage = np.zeros((len(latbins)-1,len(lonbins)-1))
        global_depth = np.zeros((len(latbins)-1,len(lonbins)-1))

        #Definition of lists to store data depending on the file
        latitudes_un = []
        longitudes_un = []
        depths_un = []


        for file in files:

            data_csv = pd.read_csv(fileDest + file , header=0, sep=',')
            latitudes_un.append(data_csv['lat'])
            longitudes_un.append(data_csv['lon'])
            depths_un.append(data_csv['depth'])

        # Loop over the global matrix
        for ii in np.arange(len(lonbins)-1):
            for jj in np.arange(len(latbins)-1):
                
                # Sub-division of each interval in between two following lat/lon into nPt_int points -> nPt-1 refinements
                long = np.linspace(lonbins[ii],lonbins[ii+1],nPt_int)
                lati = np.linspace(latbins[jj],latbins[jj+1],nPt_int)
                loc_depths = []     #list of matrices containing depths corresponding to the cells in the refined grill
                loc_bool_occs = []      #list of matrices with 1 if there exists values in the refined cell, 0 else
                #occ_matrices = []      

                for nfile in np.arange(len(files)):     # Loop over unstructured data files 
                    
                    # Local refined matrices (nRef x nRef)
                    
                    loc_sum = np.zeros((nRef, nRef))       #matrix with the sum of the values that are in each refined cell         
                    loc_depth = np.zeros((nRef, nRef))      #matrix with the average of the values that are in each refined cell 
                    loc_occ = np.zeros((nRef, nRef))        #matrix with the count of values in each refined cell
                    loc_mask_occ = np.zeros((nRef, nRef))       #matrix with 1 if there exists at least one value in the refined cell, 0 else

                    lats_un = latitudes_un[nfile]
                    lons_un = longitudes_un[nfile]
                    deps_un = depths_un[nfile]
                    
                    # Contruction of two histograms: occurrence and summation with unstructured values into structured bins
                    loc_occ, lonedges, latedges = np.histogram2d(lons_un, lats_un, bins = [long,lati])
                    loc_sum, lonedgesw, latedgesw = np.histogram2d(lons_un, lats_un, bins = [long,lati], weights=deps_un)
                    loc_occ = loc_occ.T
                    loc_sum = loc_sum.T

                    # Computation of the mean                
                    loc_depth[loc_occ != 0] = loc_sum[loc_occ != 0]/loc_occ[loc_occ != 0]
                    loc_depths.append(loc_depth)

                    # Counting occurrences
                    loc_mask_occ[loc_depth < 0]= True 
                    #loc_mask_occ[loc_occ < 1] = False
                    loc_bool_occs.append(loc_mask_occ)

                # Computation of the percentages
                loc_mask = loc_bool_occs[0]
                loc_depth = loc_depths[0]
                present_in_file = 0

                for nloc in np.arange(1,len(files)):

                    loc_mask = np.logical_or(loc_mask, loc_bool_occs[nloc])
                    if loc_bool_occs[0].sum()>0:
                        present_in_file = present_in_file + 1

                    if loc_bool_occs[nloc].sum()>0:
                        present_in_file = present_in_file + 1

                    loc_depth = np.ma.array((loc_depth, loc_depths[nloc])).sum(axis=0)

                global_percentage[jj][ii] = loc_mask.sum()/(nRef)**2
                
                if present_in_file !=0:
                    global_depth[jj][ii] = (np.nansum(loc_depth)/(nRef)**2) / present_in_file 

                else:
                    global_depth[jj][ii] = nan

                print('global percentage equal to ', global_percentage[jj][ii])
                print('global depth equal to ', global_depth[jj][ii])


        avg_depth = xr.DataArray(global_depth, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        avg_depth.loc[:,avg_depth.longitude[0]] = 0.0
        avg_depth.to_netcdf(fileDest + f'Averaged_depth_{nRef}ref.nc')
        
        perc = 50
        
        filtered_depth = xr.where(global_percentage >= perc/100, global_depth, nan)
        filtered_depth = xr.DataArray(filtered_depth, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        #saved_depth[global_percentage < perc/100] = nan
        #filtered_depth = xr.DataArray(saved_depth, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        #saved_dep_np = filtered_depth.to_numpy()
        plt.figure()
        plt.pcolormesh(filtered_depth, vmax = -0.001)
        plt.colorbar()
        plt.savefig(fileDest + f'Filtered_depths_{perc}perc_{nRef}ref', bbox_inches = 'tight', pad_inches = 0)
        plt.close() 
        filtered_depth.to_netcdf(fileDest + f'Filtered_depths_{perc}perc_{nRef}ref.nc')
        

        occ_percentage = xr.DataArray(global_percentage, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        occ_percentage_np = occ_percentage.to_numpy()
        filtered_percentage = np.ma.masked_where(occ_percentage_np < 0.4, occ_percentage_np)
        np.save(f'{fileDest}'+f'Total_filtered_percentage_{nRef}ref.npy', filtered_percentage.mask)

        plt.figure()
        plt.pcolormesh(filtered_percentage)
        plt.colorbar()

        plt.savefig(fileDest + f'Water_coverage_{perc}perc_{nRef}ref', bbox_inches = 'tight', pad_inches = 0)
        plt.close() 

        occ_percentage.to_netcdf(fileDest + f'Water_coverage_{perc}perc_{nRef}ref.nc')
        perc_dataframe = occ_percentage.to_dataframe(name = f'Water_coverage_{perc}perc')
        perc_dataframe.to_csv(fileDest + f'Water_coverage_{perc}perc_{nRef}ref.csv')

    else:
        
        avg_depth = xr.open_dataarray(f'{fileDest}' + f'Averaged_depth_{nRef}ref.nc')
        filtered_depth= xr.open_dataarray(fileDest + f'Filtered_depths_{perc}perc_{nRef}ref.nc')

    return avg_depth, filtered_depth


#####################
#       MAIN        #
#####################


fileLoc = '/g100/home/userexternal/gbuccino/Venice_Lagoon_bathymetry/'
fileDest = '/g100_scratch/userexternal/gbuccino/lagoon_analysis/Data_preparation/' 

dataset = xr.open_dataset(fileLoc + 'bathy_ADRI_CADEAU_NS.nc')
original_bathymetry = dataset['depth'] 
#original_bathymetry= xr.where(original_bathymetry == 0, nan, original_bathymetry)

plt.figure()
plt.pcolormesh(original_bathymetry, vmax = -0.001)
plt.colorbar()
plt.savefig(fileDest + f'Original_bathymetry', bbox_inches = 'tight', pad_inches = 0)
plt.close() 

NAS_bathymetry = dataset['depth'] 
lon = dataset['longitude']
lat = dataset['latitude']

plt.figure()
plt.pcolormesh(NAS_bathymetry, vmax = -0.001)
plt.colorbar()
plt.savefig(fileDest + f'NAS_Original_bathymetry', bbox_inches = 'tight', pad_inches = 0)
plt.close() 
# bathymetry.plot(vmin = -1.0, vmax =0.0)
# plt.savefig("bathy_ADRI_CADEAU_NS")

bathy_lagoon = original_bathymetry.sel(longitude = slice(12.22265625, 12.691406250), latitude= slice(45.12109375,45.589843750))
bathy_lagoon = bathy_lagoon* 100
lon_lagoon = lon.sel(longitude = slice(12.22265625, 12.691406250)).values
lat_lagoon = lat.sel(latitude= slice(45.12109375,45.589843750)).values

original_files = ['Bathy_2003CORILA', 'Bathy_2013_coarsed'] #

nPt = [3,4,5,6,7,8] #minore

converted_file_list = []
converted_namefile_list = []

for file in original_files:

    #modificare passando lista di file in cui il 2013 se ha dati sovrascrive i quadratini raffinati dell'istogramma del 2003, mentre dove non ho dati 2013 tengo 2003

    init_ref = "3004"
    fin_ref = "4326"

    unstructured_csv = pd.read_csv(f'{fileLoc}'+ file + '.csv', header=None,sep=',')
    new_name = change_ref_system(unstructured_csv, fileDest, file, init_ref, fin_ref)
    converted_file_list.append(new_name + '.csv')
    converted_namefile_list.append(new_name)

for nPt_int in nPt:

    perc = 50
    bathy_lag, filtered_depth = interp_on_structured(lon_lagoon, lat_lagoon, converted_file_list, converted_namefile_list, nPt_int, fileDest, perc)

    plt.figure()
    plt.pcolormesh(filtered_depth, vmax = -0.001)
    plt.colorbar()
    plt.savefig(fileDest + f'Filtered_depths_{perc}perc_OUT_FUNC_{nPt_int}', bbox_inches = 'tight', pad_inches = 0)
    plt.close() 
    #new_masked_bathy = xr.open_dataset(fileDest + f'Bathymetry_masked_{perc}perc_{nPt_int-1}ref.nc')
    #new_masked_bathymetry = np.ma.masked_where(filtered_depth > - 0.2, filtered_depth)
    
    # np.save(f'{fileDest}'+f'Averaged_depths_{nPt_int-1}ref.npy', new_masked_bathy.data)
    # np.save(f'{fileDest}'+f'Mask_depths_{nPt_int-1}ref.npy', new_masked_bathy.mask)

        
    # fig = plt.figure()
    # plt.pcolormesh(new_masked_bathy)
    # plt.colorbar()
    # plt.savefig(fileDest + f"MITgcm_{nPt_int-1}ref_masked_avarage_bathymetry", bbox_inches = 'tight', pad_inches = 0)
    # plt.close()


    # new_masked_bathymetry.plot(vmax =-0.001)

    # plt.savefig(fileDest + f"Lagoon_avg_bathymetry_{nPt_int-1}ref", bbox_inches = 'tight', pad_inches = 0)
    # plt.close()

    plt.figure()
    plt.pcolormesh(original_bathymetry, vmax = -0.001)
    plt.colorbar()
    plt.savefig(fileDest + f'bathy-PRE_outSUM_{perc}perc__{nPt_int}', bbox_inches = 'tight', pad_inches = 0)
    plt.close() 

    with xr.set_options(arithmetic_join="outer"):
        bathymetry = filtered_depth + original_bathymetry

    plt.figure()
    plt.pcolormesh(bathymetry, vmax = -0.001)
    plt.colorbar()
    plt.savefig(fileDest + f'BATHYMETRY_{perc}perc_OUT_FUNC_{nPt_int}', bbox_inches = 'tight', pad_inches = 0)
    plt.close() 

    bathymetry= xr.where(bathymetry != bathymetry, NAS_bathymetry, bathymetry)
    #bathymetry= xr.where(bathymetry > -0.2, nan, NAS_bathymetry)
    
    bathymetry.loc[:,bathymetry.longitude[0]] = 0.0

    bathymetry.loc[bathymetry.latitude[219],bathymetry.longitude[0]:bathymetry.longitude[6]] = 0.0
    bathymetry.loc[bathymetry.latitude[220],bathymetry.longitude[0]:bathymetry.longitude[8],] = 0.0
    bathymetry.loc[bathymetry.latitude[257]:bathymetry.latitude[263],bathymetry.longitude[44]] = 0.0
    bathymetry.loc[bathymetry.latitude[263],bathymetry.longitude[45]] = 0.0
    bathymetry.loc[bathymetry.latitude[257]:bathymetry.latitude[263],bathymetry.longitude[46]] = 0.0

    #bathymetry = np.ma.masked_where(bathymetry > - 0.2, bathymetry)

    # x = [0:6], y = 219
    # x = [0:8], y = 220
    # x = 44, y = [257, 264]
    # x = 45, y = 264
    # x = 46, y = [257, 264]

    bathymetry.to_netcdf(fileDest + f'NAS_lagoon_bathymetry_{nPt_int-1}ref_{perc}perc.nc')
    bathymetry.values.astype('f4').tofile(fileDest + f'NAS_lagoon_bathymetry_{nPt_int-1}ref_{perc}perc.bin')

    fig = plt.figure()
    bathymetry.plot(vmin = -1.5, vmax =-0.001)

    #plt.gcf().set_size_inches(40, 24)
    plt.savefig(fileDest + f"NAS_lagoon_bathymetry_{nPt_int-1}ref_{perc}perc", bbox_inches = 'tight', pad_inches = 0)

    plt.close()

print("ciao")
