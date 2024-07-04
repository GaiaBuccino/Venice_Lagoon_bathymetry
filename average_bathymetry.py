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


""" This script performs the addition of a bathymetry described in unstructured, external csv files to an existing structured bathymetry.
The procedure to combine the bathymetreis consist of averaging the external data on the structured grid and to post-process the result through a direct observation

Test case: In the example under analysis, the procedure is performed on a 1/128 refined grid of the North Adriatic Sea, to which the bathymetry of Venice lagoon is added. 
Data related to the lagoon come from to data files containing unstructured grids: CORILA2003 and coarsed2013 (that require pre-processing operations)
For a more detailed description of the algorithm used see: https://github.com/GaiaBuccino/Venice_Lagoon_bathymetry/blob/main/README.md
""" 

def transform_coordinates_system(Xin:np.ndarray, Yin:np.ndarray, EPSG_0:str ="EPSG:4326", EPSG_F:str="EPSG:3004"):
    """Function that given a meshgrid object converts its coordinates from a reference system to another

    Args:
        Xin (np.ndarray): x-coordinates coming from a meshgrid structure expressed wrt the initial reference system
        Yin (np.ndarray): y-coordinates coming from a meshgrid structure expressed wrt the initial reference system
        EPSG_0 (str, optional): starting reference system. Defaults to "EPSG:4326".
        EPSG_F (str, optional): final reference system. Defaults to "EPSG:3004".

    Returns:
        X (np.ndarray): x-coordinates expressed wrt the initial reference system
        Y (np.ndarray): y-coordinates expressed wrt the initial reference system
    """

    import pandas as pd
    import geopandas
    from shapely.geometry import Point, Polygon
    (Xflat,Yflat)=(Xin.flatten(),Yin.flatten())
    input_flatgrid=np.vstack((Xflat,Yflat)).transpose()
    if(EPSG_0=="EPSG:4326"):
        df = pd.DataFrame(input_flatgrid, columns = ['Lon','Lat'])
    else:
        df = pd.DataFrame(input_flatgrid, columns = ['Xorig','Yorig'])
    df['geometry'] = df.apply(lambda row: Point(row[df.columns[0]], row[df.columns[1]]), axis=1)
    geop_df  = geopandas.GeoDataFrame(df)
    geop_df.crs=EPSG_0 # input coordinates
    geop_df=geop_df.to_crs(EPSG_F)  # output coordinates\
    X=geop_df.geometry.x.to_numpy(dtype=float).reshape(np.shape(Xin))
    Y=geop_df.geometry.y.to_numpy(dtype=float).reshape(np.shape(Yin))
    print('old XY (',EPSG_0,') had min,average,max: X(',np.min(Xin),np.average(Xin),np.max(Xin),') Y(',np.min(Yin),np.average(Yin),np.max(Yin),') and shape',np.shape(Xin),np.shape(Yin))
    print('new XY (',EPSG_F,') had min,average,max: X(',np.min(X),np.average(X),np.max(X),') Y(',np.min(Y),np.average(Y),np.max(Y),') and shape',np.shape(X),np.shape(Y))
    return X,Y 

def geom_creation(coord1: np.ndarray, coord2:np.ndarray, structured:bool = False): 

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

                pt = Point(ii,jj)    # Coordinates of points given as lat, lon
                points.append(pt)

        return points

    else:
        pandaspoints = pd.DataFrame({'lon':coord1, 'lat':coord2})
        pandaspoints['coords'] = list(zip(pandaspoints['lon'],pandaspoints['lat']))
        pandaspoints['coords'] = pandaspoints['coords'].apply(Point)
        return pandaspoints['coords']

def change_reference(pd_bathy: pd.DataFrame, fileDest: str , fileName:str, initial_ref: str, final_ref: str):
    """Function that converts a pandas DataFrame into a GeoDataFrame structure (the geometry variable is added) and change reference system 

    Args:
        pd_bathy (pd.DataFrame): unstructured csv file with longitude, latitude, depth
        fileDest (str): path representing where the produced files have to be saved
        fileName (str): name of the unstructured file (without extension)
        initial_ref (str): starting reference system with only the code (if the reference system is EPSG:3004 -> initial_ref is '3004')
        final_ref (str): final reference system with the same rule as initial_ref

    Returns:
        (str): name of the converted file -> {fileName}_{final_ref}
    """

    new_ref_path = fileDest +f'{fileName}_{final_ref}.csv'
    
    if not os.path.exists(new_ref_path):
        
        coord_start = "EPSG:"+ initial_ref
        geom = geom_creation(pd_bathy[0], pd_bathy[1])
        gpd_bathy = gpd.GeoDataFrame(pd_bathy, geometry = geom, crs = coord_start ) 
        gpd_bathy_nr = gpd_bathy.to_crs(final_ref)      # nr -> new reference
        lat_nr= np.array(gpd_bathy_nr['geometry'].y)
        lon_nr = np.array(gpd_bathy_nr['geometry'].x)

        # Construction of the DataFrames with only necessary info 
        
        if fileName == 'Bathy_2003CORILA':              # First two columns represent the coordinates (lon, lat), third column is the depth wrt IGM 1942
            
            depth_IGM = gpd_bathy_nr[2]
            
        elif fileName == 'Bathy_2013_coarsed':

            #depth_IGM = gpd_bathy_nr[2] + 0.0243        # Conversion of data taken wrt Punta Salute into data wrt IGM 1942
            depth_IGM = gpd_bathy_nr[2] - 0.243        # Conversion of data taken wrt Punta Salute into data wrt IGM 1942

        gpd_df = {"lon": lon_nr, "lat": lat_nr, "depth": depth_IGM}
        gpd_bathy_nr = pd.DataFrame(gpd_df)
        
        # Set to zero positive depths to properly compute the average

        gpd_bathy_nr.loc[gpd_bathy_nr["depth"] > 0.0, "depth"] = 0.0  
        gpd_bathy_nr.to_csv(fileDest +f'{fileName}_{final_ref}.csv', header= ['lon','lat','depth'], index=False)

    return f'{fileName}_{final_ref}' 

def cleaning(bathy: xr.DataArray, percentage: int, coverage_percentage: np.ndarray, operations:bool=False): 
    """AD HOC FUNCTION to remove critical points

    Args:
        bathy (xr.DataArray): variable containing a map of the values to be selected
        percentage (int): lower bound of the acceptable percentage of water present in each cell
        coverage_percentage (np.ndarray): variable containing a map of the percentage of water covering each cell
        operations (bool, optional): boolean variable that express the reason why each cell is eliminated. Defaults to False.
                                     

    Returns:
        [bathy, opertions] (xr.DataArray): bathy is the map containing the values of the cleaned bathymetry
                                           operations map why the operation made on each cell (0 = not touched, != 0 removed and why)
                                                Possible values and meanings:   0 PRESENT CELL: not touched
                                                                                1 REMOVED CELL: not significant
                                                                                2 REMOVED CELL: water coverage <20%
                                                                                3 REMOVED CELL: depth value given
                                                                                4 PRESENT CELL: depth value received
                                                                                5 PRESENT CELL: depth value manually changed
                                                                                6 REMOVED CELL: superposition with NA NAS_bathymetry
                                                                                7 REMOVED CELL: depth less than hFac
                                                                                8 REMOVED CELL: closure of the domain  
                                                                                9 PRESENT CELL: averaged depths between adiacent cells
    """

    hFac = -0.2
    val = nan

    operation = xr.DataArray(np.zeros((len(bathy.longitude), len(bathy.latitude))), dims = ["longitude","latitude"], coords = {"longitude": bathy.longitude, "latitude":  bathy.latitude,})

    if operations:
        
        operation = xr.where(bathy < hFac, operation, 7)
        operation = xr.where(coverage_percentage >= percentage/100, operation, 2)
        operation.loc[:,operation.longitude[0]] = 8
        operation.loc[operation.latitude[53],operation.longitude[47]:operation.longitude[52]] = 1
        operation.loc[operation.latitude[47]:operation.latitude[52],operation.longitude[44]:operation.longitude[53]] = 1
        operation.loc[operation.latitude[54],operation.longitude[47]] = 1
        operation.loc[operation.latitude[38]:operation.latitude[39],operation.longitude[22]:operation.longitude[23]] = 1
        operation.loc[operation.latitude[40]:operation.latitude[41],operation.longitude[26]] = 1
        operation.loc[operation.latitude[28],operation.longitude[14]] = 1
        operation.loc[operation.latitude[26],operation.longitude[13]] = 1
        operation.loc[operation.latitude[15]:operation.latitude[16],operation.longitude[11]] = 1
        operation.loc[operation.latitude[9],operation.longitude[5]:operation.longitude[7]] = 1
        operation.loc[operation.latitude[8],operation.longitude[4]:operation.longitude[6]] = 1

        # With redirected channel (Pellestrina)
        
        bathy.loc[bathy.latitude[15]:bathy.latitude[19],bathy.longitude[9]] = 0
        bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[10]] = 0
        bathy.loc[bathy.latitude[16]:bathy.latitude[20],bathy.longitude[8]] = 9
        bathy.loc[bathy.latitude[20],bathy.longitude[9]:bathy.latitude[10]] = 9
        bathy.loc[bathy.latitude[21]:bathy.latitude[23],bathy.longitude[10]] = 9
        bathy.loc[bathy.latitude[23]:bathy.latitude[25],bathy.longitude[11]] = 9
        bathy.loc[bathy.latitude[25]:bathy.latitude[26],bathy.longitude[12]] = 9

        # Venezia

        operation.loc[operation.latitude[41],operation.longitude[12]:operation.longitude[15]] = 1
        operation.loc[operation.latitude[39],operation.longitude[14]] = 1
        operation.loc[operation.latitude[40],operation.longitude[14]] = 1        
        operation.loc[operation.latitude[39],operation.longitude[12]:operation.longitude[13]] = 4
        operation.loc[operation.latitude[39],operation.longitude[15]:operation.longitude[17]] = 4
        operation.loc[operation.latitude[40],operation.longitude[12]:operation.longitude[13]] = 3
        operation.loc[operation.latitude[40],operation.longitude[15]:operation.longitude[17]] = 3

        # Murano

        operation.loc[operation.latitude[43],operation.longitude[16]:operation.longitude[17]] = 1

        # Sant'Erasmo

        operation.loc[operation.latitude[41],operation.longitude[23]] = 5
        operation.loc[operation.latitude[44],operation.longitude[22]] = 5
        operation.loc[operation.latitude[42],operation.longitude[22]:operation.longitude[23]] = 1
        operation.loc[operation.latitude[43],operation.longitude[23]:operation.longitude[25]] = 1
        operation.loc[operation.latitude[44],operation.longitude[25]:operation.longitude[26]] = 1
        operation.loc[operation.latitude[45],operation.longitude[26]] = 1
        operation.loc[operation.latitude[44]:operation.latitude[45],operation.longitude[27]:operation.longitude[28]] = 9
        operation.loc[operation.latitude[43],operation.longitude[26]:operation.longitude[27]] = 9
        operation.loc[operation.latitude[42],operation.longitude[25]:operation.longitude[26]] = 9
        operation.loc[operation.latitude[41],operation.longitude[24]:operation.longitude[25]] = 9
        operation.loc[operation.latitude[39]:operation.latitude[40],operation.longitude[24]:operation.longitude[25]] = 9

        # Small channels

        operation.loc[operation.latitude[31]:operation.latitude[45],operation.longitude[0]:operation.longitude[3]] = 1
        operation.loc[operation.latitude[40],operation.longitude[4]] = 1
        operation.loc[operation.latitude[22],operation.longitude[10]] = 4
        operation.loc[operation.latitude[22],operation.longitude[11]] = 3
        operation.loc[operation.latitude[21],operation.longitude[10]] = 3
        operation.loc[operation.latitude[21],operation.longitude[9]] = 4
        operation.loc[operation.latitude[44]:operation.latitude[45],operation.longitude[6]] = 1
        operation.loc[operation.latitude[45]:operation.latitude[46],operation.longitude[7]] = 1

        # N-E single cells 

        operation.loc[operation.latitude[55]:operation.latitude[57],operation.longitude[39]:operation.longitude[41]] = 1
        operation.loc[operation.latitude[50]:operation.latitude[54],operation.longitude[42]:operation.longitude[45]] = 1
        operation.loc[operation.latitude[45]:operation.latitude[47],operation.longitude[34]:operation.longitude[43]] = 1
        operation.loc[operation.latitude[43]:operation.latitude[44],operation.longitude[30]:operation.longitude[33]] = 1

        # Superposition with NAS

        operation.loc[operation.latitude[12]:operation.latitude[16],operation.longitude[10]:operation.longitude[15]] = 6
        operation.loc[operation.latitude[27],operation.longitude[15]] = 6
        operation.loc[operation.latitude[36]:operation.latitude[40],operation.longitude[28]:operation.longitude[32]] = 6
        operation.loc[operation.latitude[36]:operation.latitude[38],operation.longitude[25]:operation.longitude[27]] = 6
        operation.loc[operation.latitude[37],operation.longitude[24]] = 6
   

    bathy = xr.where(coverage_percentage >= percentage/100, bathy, nan)
    bathy = xr.where(bathy < hFac, bathy, nan) 
    bathy.loc[:,bathy.longitude[0]] = val
    bathy.loc[bathy.latitude[53],bathy.longitude[47]:bathy.longitude[52]] = val
    bathy.loc[bathy.latitude[47]:bathy.latitude[52],bathy.longitude[44]:bathy.longitude[53]] = val
    bathy.loc[bathy.latitude[54],bathy.longitude[47]] = val
    bathy.loc[bathy.latitude[38]:bathy.latitude[39],bathy.longitude[22]:bathy.longitude[23]] = val
    
    bathy.loc[bathy.latitude[40]:bathy.latitude[41],bathy.longitude[26]] = val
    bathy.loc[bathy.latitude[28],bathy.longitude[14]] = val
    bathy.loc[bathy.latitude[26],bathy.longitude[13]] = val
    bathy.loc[bathy.latitude[15]:bathy.latitude[16],bathy.longitude[11]] = val
    bathy.loc[bathy.latitude[9],bathy.longitude[5]:bathy.longitude[7]] = val
    bathy.loc[bathy.latitude[8],bathy.longitude[4]:bathy.longitude[6]] = val

    # With redirected channel (Pellestrina)

    flow_rate = bathy.loc[bathy.latitude[15]:bathy.latitude[19],bathy.longitude[9]].sum() + bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[10]].sum()
    avg_flow_rate = flow_rate/15        # 15 number of cells that receive water from Pellestrina channel

    bathy.loc[bathy.latitude[15]:bathy.latitude[19],bathy.longitude[9]] = val
    bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[10]] = val
    bathy.loc[bathy.latitude[16]:bathy.latitude[20],bathy.longitude[8]] = bathy.loc[bathy.latitude[16]:bathy.latitude[20],bathy.longitude[8]] + avg_flow_rate
    bathy.loc[bathy.latitude[20],bathy.longitude[9]:bathy.latitude[10]] = bathy.loc[bathy.latitude[20],bathy.longitude[9]:bathy.latitude[10]] + avg_flow_rate
    bathy.loc[bathy.latitude[21]:bathy.latitude[23],bathy.longitude[10]] = bathy.loc[bathy.latitude[21]:bathy.latitude[23],bathy.longitude[10]] + avg_flow_rate
    bathy.loc[bathy.latitude[23]:bathy.latitude[25],bathy.longitude[11]] = bathy.loc[bathy.latitude[23]:bathy.latitude[25],bathy.longitude[11]] +avg_flow_rate
    bathy.loc[bathy.latitude[25]:bathy.latitude[26],bathy.longitude[12]] = bathy.loc[bathy.latitude[25]:bathy.latitude[26],bathy.longitude[12]] +avg_flow_rate
    # bathy.loc[bathy.latitude[16],bathy.longitude[9]] = val
    # bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[8]] = bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[8]] + bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[9]]+ bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[10]]
    # bathy.loc[bathy.latitude[17]:bathy.latitude[19],bathy.longitude[9]:bathy.longitude[10]] = val
    # bathy.loc[bathy.latitude[20],bathy.longitude[9]] = bathy.loc[bathy.latitude[20],bathy.longitude[9]] + bathy.loc[bathy.latitude[20],bathy.longitude[10]]
    # bathy.loc[bathy.latitude[20],bathy.longitude[10]] = val

    # Venezia

    bathy.loc[bathy.latitude[41],bathy.longitude[12]:bathy.longitude[15]] = val
    bathy.loc[bathy.latitude[39],bathy.longitude[12]:bathy.longitude[13]] = bathy.loc[bathy.latitude[39],bathy.longitude[12]:bathy.longitude[13]] + bathy.loc[bathy.latitude[40],bathy.longitude[12]:bathy.longitude[13]]
    bathy.loc[bathy.latitude[39],bathy.longitude[15]:bathy.longitude[17]] = bathy.loc[bathy.latitude[39],bathy.longitude[15]:bathy.longitude[17]] + bathy.loc[bathy.latitude[40],bathy.longitude[15]:bathy.longitude[17]]
    bathy.loc[bathy.latitude[40],bathy.longitude[12]:bathy.longitude[17]] = val

    # Murano

    bathy.loc[bathy.latitude[43],bathy.longitude[16]:bathy.longitude[17]] = val

    # Sant'Erasmo

    bathy.loc[bathy.latitude[41],bathy.longitude[23]] = hFac - 0.001
    bathy.loc[bathy.latitude[44],bathy.longitude[22]] = hFac - 0.001
    bathy.loc[bathy.latitude[42],bathy.longitude[22]:bathy.longitude[23]] = val
    bathy.loc[bathy.latitude[43],bathy.longitude[23]:bathy.longitude[25]] = val
    bathy.loc[bathy.latitude[44],bathy.longitude[25]:bathy.longitude[26]] = val
    bathy.loc[bathy.latitude[45],bathy.longitude[26]] = val

    bathy.loc[bathy.latitude[44]:bathy.latitude[45],bathy.longitude[27]:bathy.longitude[28]] = bathy.loc[bathy.latitude[44]:bathy.latitude[45],bathy.longitude[27]:bathy.longitude[28]].mean()
    bathy.loc[bathy.latitude[43],bathy.longitude[26]:bathy.longitude[27]] = bathy.loc[bathy.latitude[43],bathy.longitude[26]:bathy.longitude[27]].mean()
    bathy.loc[bathy.latitude[42],bathy.longitude[25]:bathy.longitude[26]] = bathy.loc[bathy.latitude[42],bathy.longitude[25]:bathy.longitude[26]].mean()
    bathy.loc[bathy.latitude[41],bathy.longitude[24]:bathy.longitude[25]] = bathy.loc[bathy.latitude[41],bathy.longitude[24]:bathy.longitude[25]].mean()
    bathy.loc[bathy.latitude[39]:bathy.latitude[40],bathy.longitude[24]:bathy.longitude[25]] = bathy.loc[bathy.latitude[39]:bathy.latitude[40],bathy.longitude[24]:bathy.longitude[25]].mean()
    # repeating average to uniform depths
    bathy.loc[bathy.latitude[39]:bathy.latitude[40],bathy.longitude[24]:bathy.longitude[25]] = bathy.loc[bathy.latitude[39]:bathy.latitude[40],bathy.longitude[24]:bathy.longitude[25]].mean()

    # Small channels

    bathy.loc[bathy.latitude[31]:bathy.latitude[45],bathy.longitude[0]:bathy.longitude[3]] = val
    bathy.loc[bathy.latitude[40],bathy.longitude[4]] = val
    bathy.loc[bathy.latitude[22],bathy.longitude[10]] = bathy.loc[bathy.latitude[22],bathy.longitude[10]] + bathy.loc[bathy.latitude[22],bathy.longitude[11]]
    bathy.loc[bathy.latitude[22],bathy.longitude[11]] = val
    bathy.loc[bathy.latitude[21],bathy.longitude[9]] = bathy.loc[bathy.latitude[21],bathy.longitude[9]] + bathy.loc[bathy.latitude[21],bathy.longitude[10]]
    bathy.loc[bathy.latitude[44]:bathy.latitude[45],bathy.longitude[6]] = val
    bathy.loc[bathy.latitude[45]:bathy.latitude[46],bathy.longitude[7]] = val

    # N-E single cells 

    bathy.loc[bathy.latitude[55]:bathy.latitude[57],bathy.longitude[39]:bathy.longitude[41]] = val
    bathy.loc[bathy.latitude[50]:bathy.latitude[54],bathy.longitude[42]:bathy.longitude[45]] = val
    bathy.loc[bathy.latitude[45]:bathy.latitude[47],bathy.longitude[34]:bathy.longitude[43]] = val
    bathy.loc[bathy.latitude[43]:bathy.latitude[44],bathy.longitude[30]:bathy.longitude[33]] = val

    # Superposition with NAS

    bathy.loc[bathy.latitude[12]:bathy.latitude[16],bathy.longitude[10]:bathy.longitude[15]] = val
    bathy.loc[bathy.latitude[27],bathy.longitude[15]] = val
    bathy.loc[bathy.latitude[36]:bathy.latitude[40],bathy.longitude[28]:bathy.longitude[32]] = val
    bathy.loc[bathy.latitude[36]:bathy.latitude[38],bathy.longitude[25]:bathy.longitude[27]] = val
    bathy.loc[bathy.latitude[37],bathy.longitude[24]] = val

    bathy.plot(vmin = -1.5, vmax =-0.001)

    plt.savefig(fileDest + "Cleaned_average_bathymetry", bbox_inches = 'tight', pad_inches = 0)
    bathyyy = bathy.to_dataframe(name= 'depth')
    bathyy = bathyyy.reset_index()
    bathyy.to_csv(fileDest + 'Cleaned_average_bathymetry.csv')
    plt.close()
    
    return bathy, operation

def average_on_structured(lon_st:np.ndarray, lat_st:np.ndarray, files: List[pd.DataFrame], pt: int, fileDest: str, min_percentage: int):
    """Given a list of DataFrames with unstructured data (.csv), the average of their values is computed based on 
       their position wrt the initial structured grid in order to obtain an updated grid with more accurate values of depth

    Args:
        lon_st (np.ndarray): array with the longitude data of the structured grid on which the averaging has to be performed
        lat_st (np.ndarray): array with the latitude data of the structured grid on which the averaging has to be performed
        files (List[str]): list of files with the unstructured data that have to be averaged on the structured grid
        pt (int): number of points on each cell to perform the refinement of the grid (RM: nPt points means nPt-1 refinements)
        fileDest (str): string with the destination folder where the produced files are saved

    Returns:
        xr.DataArray: file with the grid values (lon, lat) and the averaged depths
    """

    nRef = pt-1

    cleaned_depths = fileDest + f'Averaged_lagoon_bathymetry.nc'

    if not (os.path.exists(cleaned_depths)): 

        # Enlarge the grid to have values in the center of the cells instead of on their vertices

        dlon = lon_st[1] - lon_st[0]
        dlat = lat_st[1] - lat_st[0]
        lonb = lon_st - 0.5*dlon
        latb = lat_st - 0.5*dlat

        # Computing sub-cell area (m^2)

        grid = np.meshgrid(lonb, latb, indexing='ij')
        x,y = transform_coordinates_system(grid[0],grid[1],EPSG_0="EPSG:4326",EPSG_F="EPSG:3004")
        dx = np.average((x[-1]- x[0]) / (len(x)-1))
        dy = np.average((y[:,-1]- y[:,0]) / (len(y[0])-1))

        subcell_area = dx*dy / nRef**2

        # Definition of an enlarged grid

        lonbins = np.zeros((len(lonb)+1))
        lonbins[0:len(lonb)] = lonb
        lonbins[-1] = np.max(lon_st) + 0.5*dlon
        latbins = np.zeros((len(latb)+1))
        latbins[0:len(latb)] = latb
        latbins[-1] = np.max(lat_st) + 0.5*dlat

        # Definition of the matrices required for the computation
        
        global_percentage = np.zeros((len(lonbins)-1,len(latbins)-1))
        global_depth = np.zeros((len(lonbins)-1,len(latbins)-1))

        # Definition of lists to store data depending on the file

        latitudes_un = []
        longitudes_un = []
        depths_un = []

        for file in files:

            data_csv = pd.read_csv(fileDest + file , header=0, sep=',')
            latitudes_un.append(data_csv['lat'])
            longitudes_un.append(data_csv['lon'])
            depths_un.append(data_csv['depth'])
            
        lon_coord = []
        lat_coord = []

        # Loop over the global matrix

        for jj in np.arange(len(latbins)-1):
            for ii in np.arange(len(lonbins)-1):
                
                # Sub-division of each interval in between two following lat/lon into pt points -> nPt-1 refinements
                
                lon_coord.append(ii)
                lat_coord.append(jj)

                # RM: ref_ = refined grid
                ref_lon= np.linspace(lonbins[ii],lonbins[ii+1],pt)
                ref_lat = np.linspace(latbins[jj],latbins[jj+1],pt)
                ref_depths_files = []                                   # List of matrices containing depths corresponding to the cells in the refined grid
                ref_presence_files = []                                 # List of matrices with 1 if there exists values in the refined cell, 0 else
                ref_occ_files = []                                      # List of matrices with 1 if there exists values in the refined cell, 0 else
                   

                for nfile in np.arange(len(files)):                     # Loop over unstructured data files 
                    
                    # Local refined (ref) matrices (nRef x nRef)
                    
                    ref_sum = np.zeros((nRef, nRef))                    # Matrix with the sum of the values that are in each refined cell         
                    ref_depth = np.zeros((nRef, nRef))                  # Matrix with the average of the values that are in each refined cell 
                    ref_occ = np.zeros((nRef, nRef))                    # Matrix with the count of values in each refined cell
                    ref_presence = np.zeros((nRef, nRef))               # Matrix with 1 if there exists at least one value in the refined cell, 0 else
                    channel_weight = np.zeros((nRef, nRef))
                    shallow_weight = np.zeros((nRef, nRef))

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
                    ref_depths_files.append(ref_depth)
                    ref_occ_files.append(ref_occ)

                    # Counting presence (construction of a boolean matrix to count filled cells)

                    ref_presence[ref_depth < 0]= True 
                    ref_presence_files.append(ref_presence)

                ### WEIGHING ###
                
                present_in_file = 0

                for nloc in np.arange(0,len(files)):
                    if ref_presence_files[nloc].sum() > 0:
                        present_in_file = present_in_file + 1
                    if(nloc == 0):
                        ref_presence_init = ref_presence_files[0]
                        ref_depth_init = ref_depths_files[0]
                        shallow_depth = ref_depths_files[0]

                    else:
                        channel_weight = (ref_occ_files[nloc] * 100 / subcell_area)    # 100 m^2 is the area "covered" by every channel observation 
                        shallow_weight = 1 - channel_weight
                        ref_presence_init = np.logical_or(ref_presence_init, ref_presence_files[nloc])
                        ref_depth_init = ref_depths_files[nloc]*channel_weight + shallow_depth*shallow_weight      


                if present_in_file !=0:                                                      
                    global_depth[jj][ii] = np.sum(ref_depth_init)/(nRef)**2                         
                else:                                                                         
                    global_depth[jj][ii] = nan                                                  

                global_percentage[jj][ii] = ref_presence_init.sum()/(nRef)**2                                        
                print('global percentage equal to ', global_percentage[jj][ii])    
                print('global depth equal to ', global_depth[jj][ii])              

           
        averaged_depth = xr.DataArray(global_depth, dims = ["latitude","longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        
        ### CLEANING ###

        averaged_depth, operations = cleaning(averaged_depth, min_percentage, global_percentage, operations=True)

        df_averaged = averaged_depth.to_dataframe(name = 'depth')
        df_avg = df_averaged.reset_index()
        df_avg.insert(3, "cell_lon", lon_coord)
        df_avg.insert(4, "cell_lat", lat_coord)
        df_avg.to_csv(fileDest + 'Averaged_lagoon_bathymetry.csv')

        plt.figure()
        plt.pcolormesh(averaged_depth)
        plt.colorbar()
        plt.savefig(fileDest + f'Averaged_lagoon_bathymetry', bbox_inches = 'tight', pad_inches = 0)
        plt.close() 
        averaged_depth.to_netcdf(fileDest + f'Averaged_lagoon_bathymetry.nc')
        
        df_operations = operations.to_dataframe(name = 'operation')
        df_op = df_operations.reset_index()
        df_op.insert(3, "cell_lon", lon_coord)
        df_op.insert(4, "cell_lat", lat_coord)
        df_op.to_csv(fileDest + 'Cleaning_operations_map.csv')

        occurrence_percentage = xr.DataArray(global_percentage, dims = ["latitude", "longitude"], coords = {"latitude":  lat_st, "longitude": lon_st})
        occurrence_percentage.to_netcdf(fileDest + f'Water_coverage_percentage.nc')
        percentage_dataframe = occurrence_percentage.to_dataframe(name = 'depth')
        perc_df = percentage_dataframe.reset_index()
        perc_df.insert(3, "cell_lon", lon_coord)
        perc_df.insert(4, "cell_lat", lat_coord)
        perc_df.to_csv(fileDest + f'Water_coverage_percentage.csv')

    else:
        
        averaged_depth= xr.open_dataarray(fileDest + f'Averaged_lagoon_bathymetry.nc')

    return  averaged_depth



#####################
#       MAIN        #
#####################


fileLoc = './Venice_Lagoon_bathymetry/'                                   # From where files are taken
fileDest = fileLoc + 'output/'                       # Where files are saved

if not (os.path.exists(fileDest)):
    os.mkdir(fileDest)

NA_bathy_path = fileLoc + 'bathy_ADRI_CADEAU_NS.nc'     # netCDF file with: longitude, latitude, depth

dataset = xr.open_dataset(NA_bathy_path)
original_bathy = dataset['depth'] 

# Plot of the original bathymetry 

plt.figure()
plt.pcolormesh(original_bathy, vmax = -0.001)
plt.colorbar()
plt.savefig(fileDest + f'original_bathymetry', bbox_inches = 'tight', pad_inches = 0)
plt.close() 

NAS_bathymetry = dataset['depth'] 
lon = dataset['longitude']
lat = dataset['latitude']

lagoon_bathy = original_bathy.sel(longitude = slice(12.22265625, 12.691406250), latitude= slice(45.12109375,45.589843750))

lagoon_bathy = lagoon_bathy * 0.0
lon_lagoon = lon.sel(longitude = slice(12.22265625, 12.691406250)).values
lat_lagoon = lat.sel(latitude= slice(45.12109375,45.589843750)).values

original_files = ['Bathy_2003CORILA', 'Bathy_2013_coarsed'] 	# Files with unstructured data to be averaged on the grid

nPt = [8]                                                       # The largest value for which each refined cell certainly contains points is 7 (RM: nRefinement = nInternalPoints - 1)

nr_files = []


### PRE-PROCESSING ###

for file in original_files:

    initial_ref= "3004"
    final_ref = "4326"

    unstructured_csv = pd.read_csv(f'{fileLoc}'+ file + '.csv', header=None,sep=',')
    new_name = change_reference(unstructured_csv, fileDest, file, initial_ref, final_ref)
    nr_files.append(new_name + '.csv')


for pt in nPt:

    minimal_percentage = 20

    ### AVERAGING ###

    filtered_depth = average_on_structured(lon_lagoon, lat_lagoon, nr_files, pt, fileDest, minimal_percentage)  

    ### ADDING ###

    with xr.set_options(arithmetic_join="outer"):
        bathymetry = filtered_depth + original_bathy

    bathymetry= xr.where(bathymetry != bathymetry, NAS_bathymetry, bathymetry)
    fig = plt.figure()
    bathymetry.plot(vmin = -1.5, vmax =-0.001)

    plt.gcf().set_size_inches(20, 12)                       # To preserve the dimensions of the image without stretching it
    plt.savefig(fileDest + f"NAS_lagoon_bathymetry", dpi= 150, bbox_inches = 'tight', pad_inches = 0)
    plt.close()
    
    ### EXPORTING ###

    bathymetry.to_netcdf(fileDest + f'NAS_lagoon_bathymetry.nc')
    bathymetry.values.astype('f4').tofile(fileDest + f'NAS_lagoon_bathymetry.bin')

    

print("*** The average bathymetry has been computed ***")
