# Venice Lagoon Bathymetry
Algorithm that updates the structured bathymetry of the North Adriatic Sea (resolution at 1/128 degree) adding the bathymetry of the Venice lagoon constructed by averaging values coming from more recent unstructured bathymetries (CORILA2003 and coarsed2013)

## Required files
 - Bathymetry of the North Adriatic Sea (netCDF format) on which the Venice lagoon bathymetry is innested, in the following referred as __original bathymetry__ (in the example _bathy_ADRI_CADEAU_NS.nc_)
 - Files (csv format) with depth data to update the original bathymetry, in the following referred as __external files__ ( in the example _Bathy_2003CORILA.csv_ and _Bathy_2013_coarsed.csv_)

## _average_bathymetry.py_
The procedure to construct the Venice lagoon updated bathymetry consists of 5 steps: __Pre-processing__ data, __Averaging__ data recorded in the _external files_, __Weighing__ data given their occurrence on the grid-cells, __Adding__ the averaged-weighted bathymetry to the original bthymetry, __Exporting__ the results in different formats.

To start the procedure, the user has to modify the path where the necessary files are stored _fileLoc_ and the one where the produced files have to be saved _fileDest_:
```
fileLoc = ''
fileDest = '' 
```

In the following each step is described in details:

### Pre-processing
The first step fulfilled is the conversion of the system of reference of coordinates of the _external files_ in order to be comparable with the _original bathymetry_ data. To do so, all the _external files_ are converted into GeoDataFrame to exploit the built-in function that automatically convert the coordinate reference system. The initial and final reference system have to be specified changing the values of the variables _init_ref_ representing the initial reference system and _fin_ref_ that is the final one
```
init_ref = "..."
fin_ref = "..."
```

### Averaging
The first step of the algorithm used for the averaging consists of the construction of a more refined grid in which each cell of the _original bathymetry_ is sub-divided into _nRef <sup>2</sup>_ sub-cells, where _nRef_ corresponds to the number of refinements desired and its a customizable parameter. Each sub-cell represents a bin of two different histograms that is constructed for each cell. The first histogram records per each sub-cell the number of observations that fall in its area, while the second compute their sum. At this point, the average value is computed dividing the value of the sum by the number of observation per each sub-cell. Once this computation is performed, the value of a further average operation between all the sub-cells is computed and the resulting value is assigned to the original cell. The result of this operation is a matrix with the same dimension as the original grid with in each cell the averaged value of depth coming from the two combined average procedures.
