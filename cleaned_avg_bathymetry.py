# importing required libraries

from math import nan
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
from shapely.geometry import Point
import matplotlib.pyplot as plt
import geopandas as gpd
import os
from typing import List

"""_summary_
""" 

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



def change_reference(pd_bathy: pd.DataFrame, fileDest: str , fileName:str, initial_ref: str, final_ref: str):
    """Function that converts a pandas DataFrame into a GeoDataFrame structure (it adds the geometry variable and the changing of the reference system is a built-in method)

    Args:
        pd_bathy (pd.DataFrame): _description_
        fileDest (str): _description_
        fileName (str): _description_
        initial_ref (str): _description_
        final_ref (str): _description_

    Returns:
        _type_: _description_
    """

    new_ref_path = fileDest +f'{fileName}_{final_ref}.csv'
    
    if not os.path.exists(new_ref_path):
        
        coord_start = "EPSG:"+ initial_ref
        geom = geom_creation(pd_bathy[1], pd_bathy[0])
        gpd_bathy = gpd.GeoDataFrame(pd_bathy, geometry = geom, crs = coord_start ) # type: ignore

        # nr -> new reference

        gpd_bathy_nr = gpd_bathy.to_crs(final_ref)
        lat_nr= np.array(gpd_bathy_nr['geometry'].y)
        lon_nr = np.array(gpd_bathy_nr['geometry'].x)

        # Construction of the DataFrames with only necessary info 
        
        if fileName == 'Bathy_2003CORILA':      # first two columns represent the coordinates (lon, lat), third column is the depth wrt IGM 1942
            
            depth_IGM = gpd_bathy_nr[2]
            # gpd_bathy_nr = gpd_bathy_nr.drop(gpd_bathy_nr.columns[3], axis =1)
            # gpd_bathy_nr = gpd_bathy_nr.drop(gpd_bathy_nr.columns[1], axis =1)
            # gpd_bathy_nr = gpd_bathy_nr.drop(gpd_bathy_nr.columns[0], axis =1)
            # gpd_bathy_nr.insert(0, "lon", lon_nr)
            # gpd_bathy_nr.insert(1, "lat", lat_nr)
            # gpd_bathy_nr.rename(columns={0: "lon", 1: "lat", 2:"depth", "geometry": "geometry"}, inplace = True)

            # Set to zero positive depths to properly compute the average
            # gpd_bathy_nr.loc[gpd_bathy_nr["depth"] > 0.0, "depth"] = 0.0
            
            # gpd_bathy_nr.to_csv(fileDest +f'{fileName}_{final_ref}.csv', header= ['lon','lat','depth', 'geometry'], index=False)

        elif fileName == 'Bathy_2013_coarsed':

            depth_IGM = gpd_bathy_nr[2] + 0.0243        #conversion of data taken wrt Punta Salute into data wrt IGM 1942

        gpd_df = {"lon": lon_nr, "lat": lat_nr, "depth": depth_IGM}
        gpd_bathy_nr = pd.DataFrame(gpd_df)
            # gpd_bathy_nr = gpd_bathy_nr.drop(gpd_bathy_nr.columns[1], axis =1)
            # gpd_bathy_nr = gpd_bathy_nr.drop(gpd_bathy_nr.columns[0], axis =1)
            # gpd_bathy_nr.insert(0, "lon", lon_nr)
            # gpd_bathy_nr.insert(1, "lat", lat_nr)
            # gpd_bathy_nr.rename(columns={0: "lon", 1: "lat", 2:"depth", "geometry": "geometry"}, inplace = True)

            # Set to zero positive depths to properly compute the average
        gpd_bathy_nr.loc[gpd_bathy_nr["depth"] > 0.0, "depth"] = 0.0  
        gpd_bathy_nr.to_csv(fileDest +f'{fileName}_{final_ref}.csv', header= ['lon','lat','depth'], index=False)

    return f'{fileName}_{final_ref}'    

def cleaning(bathy: xr.DataArray, percentage: int, coverage_percentage: np.ndarray, channel:bool = True): 
    """DA AUTOMATIZZARE ->  passare lista di intervalli?
    Args:
        bathy (xr.DataArray): _description_
        percentage (int): _description_
        coverage_percentage (np.ndarray): _description_
        channel (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """  

    hFac = -0.2

    bathy = xr.where(coverage_percentage >= percentage/100, bathy, nan)
    bathy = xr.where(bathy < hFac, bathy, nan)

    val = nan

    bathy.loc[:,bathy.longitude[0]] = val
    bathy.loc[bathy.latitude[53],bathy.longitude[47]:bathy.longitude[52]] = val
    bathy.loc[bathy.latitude[47]:bathy.latitude[52],bathy.longitude[44]:bathy.longitude[53]] = val
    bathy.loc[bathy.latitude[54],bathy.longitude[47]] = val
    bathy.loc[bathy.latitude[38]:bathy.latitude[39],bathy.longitude[22]:bathy.longitude[23]] = val
    bathy.loc[bathy.latitude[40]:bathy.latitude[41],bathy.longitude[26]] = val
    bathy.loc[bathy.latitude[28],bathy.longitude[14]] = val
    bathy.loc[bathy.latitude[26],bathy.longitude[13]] = val
    bathy.loc[bathy.latitude[15]:bathy.latitude[16],bathy.longitude[11]] = val
    bathy.loc[bathy.latitude[9],bathy.longitude[6]:bathy.longitude[7]] = val

    # without border channel (Pellestrina)
    if not channel:
        bathy.loc[bathy.latitude[15]:bathy.latitude[19],bathy.longitude[9]] = val
        bathy.loc[bathy.latitude[15]:bathy.latitude[22],bathy.longitude[10]] = val
        bathy.loc[bathy.latitude[22]:bathy.latitude[25],bathy.longitude[11]] = val
        bathy.loc[bathy.latitude[25]:bathy.latitude[26],bathy.longitude[12]] = val
        bathy.loc[bathy.latitude[26],bathy.longitude[13]] = val

    # with redirected channel (Pellestrina)
    else: 
        bathy.loc[bathy.latitude[15],bathy.longitude[9]] = val
        bathy.loc[bathy.latitude[16],bathy.longitude[8]] = bathy.loc[bathy.latitude[16],bathy.longitude[8]] + bathy.loc[bathy.latitude[16],bathy.longitude[9]]
        bathy.loc[bathy.latitude[16],bathy.longitude[9]] = val
        bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[8]] = bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[8]] + bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[9]]+ bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[10]]
        bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[9]:bathy.longitude[10]] = val
        bathy.loc[bathy.latitude[20],bathy.longitude[9]] = bathy.loc[bathy.latitude[20],bathy.longitude[9]] + bathy.loc[bathy.latitude[20],bathy.longitude[10]]
        bathy.loc[bathy.latitude[20],bathy.longitude[10]] = val

    # Venezia
    bathy.loc[bathy.latitude[41],bathy.longitude[12]:bathy.longitude[15]] = val
    bathy.loc[bathy.latitude[40],bathy.longitude[12]:bathy.longitude[17]] = val

    # Murano
    bathy.loc[bathy.latitude[43],bathy.longitude[16]:bathy.longitude[17]] = val

    # Sant'Erasmo
    bathy.loc[bathy.latitude[42],bathy.longitude[22]:bathy.longitude[23]] = val
    bathy.loc[bathy.latitude[43],bathy.longitude[23]:bathy.longitude[25]] = val
    bathy.loc[bathy.latitude[44],bathy.longitude[25]:bathy.longitude[26]] = val
    bathy.loc[bathy.latitude[45],bathy.longitude[26]] = val

    # Small channels
    bathy.loc[bathy.latitude[31]:bathy.latitude[45],bathy.longitude[0]:bathy.longitude[3]] = val
    bathy.loc[bathy.latitude[40],bathy.longitude[4]] = val
    bathy.loc[bathy.latitude[41]:bathy.latitude[45],bathy.longitude[4]:bathy.longitude[5]] = val
    bathy.loc[bathy.latitude[43]:bathy.latitude[45],bathy.longitude[6]] = val
    bathy.loc[bathy.latitude[45]:bathy.latitude[46],bathy.longitude[7]] = val

    # N-E single cells 
    bathy.loc[bathy.latitude[55]:bathy.latitude[57],bathy.longitude[39]:bathy.longitude[41]] = val
    bathy.loc[bathy.latitude[50]:bathy.latitude[54],bathy.longitude[42]:bathy.longitude[45]] = val
    bathy.loc[bathy.latitude[45]:bathy.latitude[47],bathy.longitude[34]:bathy.longitude[43]] = val
    bathy.loc[bathy.latitude[43]:bathy.latitude[44],bathy.longitude[30]:bathy.longitude[33]] = val

    fig = plt.figure()
    bathy.plot(vmin = -1.5, vmax =-0.001)

    plt.savefig(fileDest + "Cleaned_bathymetry", bbox_inches = 'tight', pad_inches = 0)
    bathyyy = bathy.to_dataframe(name= 'depth')
    bathyy = bathyyy.reset_index()
    bathyy.to_csv(fileDest + 'Cleaned_bathymetry.csv')

    plt.close()

    return bathy

def average_on_structured(lon_st:np.ndarray, lat_st:np.ndarray, files: List[pd.DataFrame], pt: int, fileDest: str, perc: int):
    """Given a list of files with unstructured data (csv), thay are interpolated on a structured grid through the computation of the spatial average on a more refined 
    grid in order to have more accurate values

    Args:
        lon_st (np.ndarray): array with the longitude data of the structured grid on which the interpolation has to be performed
        lat_st (np.ndarray): array with the latitude data of the structured grid on which the interpolation has to be performed
        files (List[str]): list of the name of the file with the unstructured data tht have to be interpolated on the structured grid
        pt (int): number of points on each cell to perform the refinement of the grid (RM: nPt points means nPt-1 refinements)
        fileDest (str): string with the destination folder where the produced files are saved

    Returns:
        xr.DataArray: file with the grid values (lon, lat) and the averaged depths
    """

    nRef = pt-1

    depth_path = fileDest + f'Averaged_depth_{nRef}ref.nc'
    percentage_path = fileDest + f'Water_coverage_{perc}perc_{nRef}ref.csv'

    if not (os.path.exists(depth_path) & os.path.exists(percentage_path)):

        # Enlarge the grid to have values in the center of the cells instead of on their vertices
        dlon = lon_st[1] - lon_st[0]
        dlat = lat_st[1] - lat_st[0]
        lonb = lon_st - 0.5*dlon
        latb = lat_st - 0.5*dlat

        # Definition of an enlarged grid
        lonbins = np.zeros((len(lonb)+1))
        lonbins[0:len(lonb)] = lonb
        lonbins[-1] = np.max(lon_st) + 0.5*dlon
        latbins = np.zeros((len(latb)+1))
        latbins[0:len(latb)] = latb
        latbins[-1] = np.max(lat_st) + 0.5*dlat

        # Definition of the matrices required for the computation
        
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
            
        lonHD=np.linspace(lonbins[0],lonbins[-1],pt*(len(lonbins)-1))                                                     # CL
        latHD=np.linspace(latbins[0],latbins[-1],pt*(len(latbins)-1))                                                     # CL
        depthHD=np.zeros((len(files),nRef*(len(latbins)-1),nRef*(len(lonbins)-1)),dtype=float)                                 # CL
        occurHD=np.zeros((len(files),nRef*(len(latbins)-1),nRef*(len(lonbins)-1)),dtype=float)                                 # CL
        
        lon_coord = []
        lat_coord = []


        # Loop over the global matrix
        for ii in np.arange(len(lonbins)-1):
            for jj in np.arange(len(latbins)-1):
                
                # Sub-division of each interval in between two following lat/lon into pt points -> nPt-1 refinements
                
                lon_coord.append(ii)
                lat_coord.append(jj)

                ref_lon= np.linspace(lonbins[ii],lonbins[ii+1],pt)
                ref_lat = np.linspace(latbins[jj],latbins[jj+1],pt)

                ref_depths_files = []     #list of matrices containing depths corresponding to the cells in the refined grill
                ref_presence_files = []      #list of matrices with 1 if there exists values in the refined cell, 0 else
                   

                for nfile in np.arange(len(files)):     # Loop over unstructured data files 
                    
                    # Local refined (ref) matrices (nRef x nRef)
                    
                    ref_sum = np.zeros((nRef, nRef))       #matrix with the sum of the values that are in each refined cell         
                    ref_depth = np.zeros((nRef, nRef))      #matrix with the average of the values that are in each refined cell 
                    ref_occ = np.zeros((nRef, nRef))        #matrix with the count of values in each refined cell
                    ref_presence = np.zeros((nRef, nRef))       #matrix with 1 if there exists at least one value in the refined cell, 0 else

                    lats_un = latitudes_un[nfile]
                    lons_un = longitudes_un[nfile]
                    deps_un = depths_un[nfile]
                    
                    # Contruction of two histograms: occurrence and summation with unstructured values into structured bins
                    ref_occ, lonedges, latedges = np.histogram2d(lons_un, lats_un, bins = [ref_lon,ref_lat])
                    ref_sum, lonedgesw, latedgesw = np.histogram2d(lons_un, lats_un, bins = [ref_lon,ref_lat], weights=deps_un)
                    ref_occ = ref_occ.T
                    ref_sum = ref_sum.T

                    # Computation of the mean                
                    ref_depth[ref_occ != 0] = ref_sum[ref_occ != 0]/ref_occ[ref_occ != 0]
                    depthHD[nfile,jj*nRef:(jj+1)*nRef,ii*nRef:(ii+1)*nRef] = ref_depth[:,:]                                     # CL 
                    occurHD[nfile,jj*nRef:(jj+1)*nRef,ii*nRef:(ii+1)*nRef] = ref_occ[:,:]                                       # CL
                    ref_depths_files.append(ref_depth)

                    # Counting presence
                    ref_presence[ref_depth < 0]= True 
                    ref_presence_files.append(ref_presence)

                # Computation of the percentages
                bool_occ = ref_presence_files[0]
                bool_occ_sum = ref_presence_files[0]                                               ## CL
                ref_depth = ref_depths_files[0]
                present_in_file = 0

                for nloc in np.arange(1,len(files)):

                    bool_occ = np.logical_or(bool_occ, ref_presence_files[nloc])

                    if nloc == 1 and ref_presence_files[0].sum() > 0:                                              ## CL : added if(nloc ==1)
                        
                        present_in_file = present_in_file + 1

                    if ref_presence_files[nloc].sum() > 0:
                        
                        present_in_file = present_in_file + 1

                    ref_depth = np.ma.array((ref_depth, ref_depths_files[nloc])).sum(axis=0)


                    bool_occ_sum = bool_occ_sum + ref_presence_files[nloc]                             ## CL : alternative
                
                ref_depth = ref_depth / bool_occ_sum                                              ## CL : alternative

                if present_in_file !=0:                                                       ## CL : alternative
                  
                  global_depth[jj][ii] = np.sum(ref_depth)/(nRef)**2                         ## CL : alternative
                
                else:                                                                         ## CL : alternative
                  
                  global_depth[jj][ii] = nan                                                  ## CL : alternative

                global_percentage[jj][ii] = bool_occ.sum()/(nRef)**2                           

                
                print('global percentage equal to ', global_percentage[jj][ii])    
                print('global depth equal to ', global_depth[jj][ii])              

            print("ciao")

        percent = 20
        avg_depth = xr.DataArray(global_depth, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        avg_depth.loc[:,avg_depth.longitude[0]] = 0.0
        avg_depth_cleaned = cleaning(avg_depth, percent, global_percentage)
        avg_depth_cleaned.to_netcdf(fileDest + f'Averaged_depth_{nRef}ref.nc')
     
        filtered_depth = xr.where(global_percentage >= percent/100, global_depth, nan)
        filtered_depth = xr.DataArray(filtered_depth, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        filtered_depth = cleaning(filtered_depth, percent, global_percentage)
        
       
        plt.figure()
        plt.pcolormesh(filtered_depth, vmax = -0.001)
        plt.colorbar()
        plt.savefig(fileDest + f'Filtered_depths_{perc}perc_{nRef}ref', bbox_inches = 'tight', pad_inches = 0)
        plt.close() 
        filtered_depth.to_netcdf(fileDest + f'Filtered_depths_{perc}perc_{nRef}ref.nc')

        filt = filtered_depth.to_dataframe(name = 'depth')
        filt_df = filt.reset_index()
        filt_df.to_csv(fileDest + 'FILTERED.csv')
        
        occ_percentage = xr.DataArray(global_percentage, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        # occ_percentage_np = occ_percentage.to_numpy()
        # #filtered_percentage = xr.where(occ_percentage_np < 0.4, occ_percentage_np, )
        # occ_percentage = cleaning(occ_percentage)
        # np.save(f'{fileDest}'+f'Total_filtered_percentage_{nRef}ref.npy')#, filtered_percentage.mask)

        # plt.figure()
        # plt.pcolormesh(filtered_percentage)
        # plt.colorbar()

        # plt.savefig(fileDest + f'Water_coverage_{perc}perc_{nRef}ref', bbox_inches = 'tight', pad_inches = 0)
        # plt.close() 

        occ_percentage.to_netcdf(fileDest + f'Water_coverage_{perc}perc_{nRef}ref.nc')
        perc_dataframe = occ_percentage.to_dataframe(name = 'depth')
        perc_df = perc_dataframe.reset_index()
        perc_df.insert(3, "cell_lon", lon_coord)
        perc_df.insert(4, "cell_lat", lat_coord)
        perc_df.to_csv(fileDest + f'Water_coverage_{perc}perc_{nRef}ref.csv')

    else:
        
        avg_depth = xr.open_dataarray(f'{fileDest}' + f'Averaged_depth_{nRef}ref.nc')
        filtered_depth= xr.open_dataarray(fileDest + f'Filtered_depths_{perc}perc_{nRef}ref.nc')

    return avg_depth, filtered_depth





#####################
#       MAIN        #
#####################


fileLoc = '/g100/home/userexternal/gbuccino/Venice_Lagoon_hydrodynamics/'   #from where files are taken
fileDest = '/g100_scratch/userexternal/gbuccino/lagoon_analysis/Data_preparation/'  #where files are saved
NA_bathy_path = fileLoc + 'bathy_ADRI_CADEAU_NS.nc'     # netCDF file with: longitude, latitude, depth

dataset = xr.open_dataset(NA_bathy_path)
original_bathy = dataset['depth'] 

# Plot of the original bathymetry 

plt.figure()
plt.pcolormesh(original_bathy, vmax = -0.001)
plt.colorbar()
plt.savefig(fileDest + f'original_bathymetry', bbox_inches = 'tight', pad_inches = 0)
plt.close() 

# Make a copy of the original bathymetry

NAS_bathymetry = dataset['depth'] 
lon = dataset['longitude']
lat = dataset['latitude']

lagoon_bathy = original_bathy.sel(longitude = slice(12.22265625, 12.691406250), latitude= slice(45.12109375,45.589843750))
lagoon_bathy = lagoon_bathy * 0.0
lon_lagoon = lon.sel(longitude = slice(12.22265625, 12.691406250)).values
lat_lagoon = lat.sel(latitude= slice(45.12109375,45.589843750)).values

original_files = ['Bathy_2003CORILA', 'Bathy_2013_coarsed'] 	# files with unstructured data to be averaged on the grid

nPt = [8]#np.arange(5,8,1)     # number of internal points [5,6,7,8] || The largest value for which each refinement certainly contains points is 7 (8 internal points)

nr_files = []
#nr_namefiles = []

for file in original_files:

    #modificare passando lista di file in cui il 2013 se ha dati sovrascrive i quadratini raffinati dell'istogramma del 2003, mentre dove non ho dati 2013 tengo 2003

    initial_ref= "3004"
    final_ref = "4326"

    unstructured_csv = pd.read_csv(f'{fileLoc}'+ file + '.csv', header=None,sep=',')
    new_name = change_reference(unstructured_csv, fileDest, file, initial_ref, final_ref)
    nr_files.append(new_name + '.csv')
    #nr_namefiles.append(new_name)

for pt in nPt:

    minimal_percentage = 20
    bathy_lag, filtered_depth = average_on_structured(lon_lagoon, lat_lagoon, nr_files, pt, fileDest, minimal_percentage)

    # plt.figure()
    # plt.pcolormesh(filtered_depth, vmax = -0.001)
    # plt.colorbar()
    # plt.savefig(fileDest + f'Filtered_depths_{perc}perc_OUT_FUNC_{pt}', bbox_inches = 'tight', pad_inches = 0)
    # plt.close() 
   
    # plt.figure()
    # plt.pcolormesh(original_bathy, vmax = -0.001)
    # plt.colorbar()
    # plt.savefig(fileDest + f'bathy-PRE_outSUM_{perc}perc__{pt}', bbox_inches = 'tight', pad_inches = 0)
    # plt.close() 

    # with xr.set_options(arithmetic_join="outer"):
    #     bathymetry = filtered_depth + original_bathy

    # plt.figure()
    # plt.pcolormesh(bathymetry, vmax = -0.001)
    # plt.colorbar()
    # plt.savefig(fileDest + f'BATHYMETRY_{perc}perc_OUT_FUNC_{pt}', bbox_inches = 'tight', pad_inches = 0)
    # plt.close() 

    # bathymetry= xr.where(bathymetry != bathymetry, NAS_bathymetry, bathymetry)
    # #bathymetry= xr.where(bathymetry > -0.2, nan, NAS_bathymetry)
    
    """    
    NON SERVONO PIU'

    # bathymetry.loc[bathymetry.latitude[219],bathymetry.longitude[0]:bathymetry.longitude[6]] = 0.0
    # bathymetry.loc[bathymetry.latitude[220],bathymetry.longitude[0]:bathymetry.longitude[8],] = 0.0
    # bathymetry.loc[bathymetry.latitude[257]:bathymetry.latitude[263],bathymetry.longitude[44]] = 0.0
    # bathymetry.loc[bathymetry.latitude[263],bathymetry.longitude[45]] = 0.0
    # bathymetry.loc[bathymetry.latitude[257]:bathymetry.latitude[263],bathymetry.longitude[46]] = 0.0

    # #bathymetry = np.ma.masked_where(bathymetry > - 0.2, bathymetry)

    # # x = [0:6], y = 219
    # # x = [0:8], y = 220
    # # x = 44, y = [257, 264]
    # # x = 45, y = 264
    # # x = 46, y = [257, 264] """

    # bathymetry.to_netcdf(fileDest + f'NAS_lagoon_bathymetry_{pt-1}ref_{perc}perc.nc')
    # bathymetry.values.astype('f4').tofile(fileDest + f'NAS_lagoon_bathymetry_{pt-1}ref_{perc}perc.bin')

    # fig = plt.figure()
    # bathymetry.plot(vmin = -1.5, vmax =-0.001)

    # #plt.gcf().set_size_inches(40, 24)
    # plt.savefig(fileDest + f"NAS_lagoon_bathymetry_{pt-1}ref_{perc}perc", bbox_inches = 'tight', pad_inches = 0)

    # plt.close()

print("ciao")