# Venice Lagoon Bathymetry
Algorithm that updates the structured bathymetry of the North Adriatic Sea (resolution at 1/128 degree) adding the bathymetry of the Venice lagoon constructed by averaging values coming from more recent unstructured bathymetries (CORILA2003 and coarsed2013)

## Required files
 - Bathymetry of the North Adriatic Sea (netCDF format) on which the Venice lagoon bathymetry is innested, in the following referred as __original bathymetry__ (in this example _bathy_ADRI_CADEAU_NS.nc_)
 - Files (csv format) with depth data to update the original bathymetry, in the following referred as __external files__ (in this example _Bathy_2003CORILA.csv_ and _Bathy_2013_coarsed.csv_)

## _average_bathymetry.py_
The procedure to construct the Venice lagoon updated bathymetry consists of 5 steps: __Pre-processing__ data, __Averaging__ data recorded in the _external files_, __Weighing__ data given their occurrence on the grid-cells, __Cleaning__ of critical points (clusters of solated points, land points, boundary points), __Adding__ the averaged-weighted bathymetry to the original bthymetry, __Exporting__ the results in different formats.

To start the procedure, the user has to modify the path where the necessary files are stored _fileLoc_ and the one where the produced files have to be saved _fileDest_: (they are both defined as the current folder by default).

In the following each step is described in details:

### PRE-PROCESSING 
The first step fulfilled is the conversion of the reference system of coordinates of the _external files_ in order to be comparable with the _original bathymetry_ data. To do so, all the _external files_ are converted into GeoDataFrame to exploit the built-in function that automatically convert the coordinate reference system (to_crs()). The initial and final reference system have to be specified changing the values of the variables _init_ref_ representing the initial reference system and _fin_ref_ that is the final one

>init_ref = "3004"
>
>fin_ref = "4326"

(_SEE THE FINAL REMARK FOR DETAILS ABOUT DATA_)

### AVERAGING
The first step of the algorithm used for the averaging consists of the construction of a more refined grid in which each cell of the _original bathymetry_ is sub-divided into _nRef <sup>2</sup>_ sub-cells, where _nRef_ corresponds to the number of refinements desired and its a customizable parameter.

> nRef = 7 
> 
> This value depends on the density of samplings contained in the unstructured files, it is the largest value to certainly capture external data in each sub-cells (in area where data are present).

Each sub-cell represents a bin of two different histograms that are constructed for each cell. The first histogram records per each sub-cell the number of observations that fall in its area, while the second compute the sum of the latters. At this point, the average value is computed dividing the value of the sum by the number of observation per each sub-cell. Once this computation is performed, the value of a further average (sum of the values in each sub-cell / nRef <sup>2</sup>) operation between all the sub-cells is computed and the resulting value is assigned to the original cell. The output of this operation is a matrix with the same dimension as the original grid containing in each cell the averaged value of depth coming from the two combined average procedures.

### WEIGHING 
The weighing operation is performed to take into account the distribution of data in the decision of the caracterization of each cell (water or land), the depth values computed in the averaging step is weighed with respect to the percentage of occurrence of water into each cell. In order to do so, in each cell and for each _external file_ a matrix of boolean representing the occurreces of data considering every sub-cell is constructed (_True_ values means that at least one values falls into the cell, _False_ is assigned in the opposite case) and the union of these matrices with respect to the different _external file_ is computed. The resulting matrix is used to compute the percentage of water covering each cell. If the value computed is under a certain threshold, the cell is converted into land (depth equal to zero). In this way, all the cells that have only a small area covered by water are neglected into the reconstruction of the bathymetry (this allows to avoid enlargement of the surface covered by water through very shallow areas)

### CLEANING
The cleaning procedure consists in the visual correction of the outcome of the previous operations in order to have a result consistent with the reality.
Each cell is categorized into a categorical variable as follows:

| Syntax      | Description | Test Text     |
| :---        |    :----:   |          ---: |
| Header      | Title       | Here's this   |
| Paragraph   | Text        | And more      |

- **0**  -> not touched
>
>**1**  REMOVED CELL: not significant
>
>**2**  REMOVED CELL: water coverage < 20% (customizable value changing the value of the variable _min_percentage_)
>
>**3**  REMOVED CELL: depth value given
>
>**4** REMOVED CELL: depth value received
>
>**5** REMOVED CELL: depth value manually changed
>
>**6** REMOVED CELL: superposition with original bathymetry
>
>**7** REMOVED CELL: depth less than hFac
>
>**8** REMOVED CELL: closure of the domain

This operation is ad hoc for each case and is performed introducing 
>__operation__.loc[operation.latitude[_latitude index_],operation.longitude[_longitude index_]] = _correspondent value_
>__operation__.loc[operation.latitude[_initial_latitude index_]:operation.latitude[_final_latitude index_],operation.longitude[_initial_longitude index_]:operation.longitude[_final_longitude index_]] = _correspondent value_
>
The same modifications has to be done on the __bathy__ variable

<ins>Remark:</ins> values 3 and 4 have been introduced to guarantee the conservation of the volume of water. This is done adding the value of the depth of some particular cells to the one adiacent, based on the necessity to avoid the connection of cells or to represent land cells without neglecting important water channels

### ADDING
Once the final bathymetry is ready, it is added to the _original bathymetry_ of the North Adriatic Sea given at the beginning through the simple superposition of the values present in the two different files, the _original bathymetry_ and the _lagoon bathymetry_

### EXPORTING
The produced files are exported in the .nc and .bin format

### Final Remarks
<ins>About Data</ins>

- _original structured bathymetry:_ data are collected with respect to the Geographic Military Institute zero, called Genova 1942 (IGM 1942)
- _Bathy_CORILA2003.csv:_ contains data collected with respect to the Geographic Military Institute zero in the third column (index 2 in python) and with respect to the zero of Punta Salute in the last column
- _Bathy_2013_coarsed.csv:_ data are referred to the zero of Punta Salute, so they need to be rescaled to be referred to the zero of IGM 1942. This is done by adding 0.243 to each value of depth 

_see Bathymetry_scaling.pdf for details and visual explanation_
