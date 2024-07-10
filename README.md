# Venice Lagoon Bathymetry
Algorithm that updates the bathymetry of the structured mesh of the North Adriatic Sea adding the bathymetry of the Venice lagoon constructed by averaging values coming from more recent unstructured data :

- **CORILA2003** : dataset collected using on a “Multibeam” bathymetric acquisition system, from the end of 1999 to spring 2003, covering the whole lagoon including the most shallow areas, and not only the navigable channels. This bathymetry dataset reports depths with respect to the *zero IGM 1942*. In this file, coordinates are expressed in meters


- **coarsed2013** : high-resolution ASCII:ESRI gridded bathymetric data (0.5 m DTM) [Madricardo et al, 2016](http://dx.doi.org/10.1594/IEDA/323605), [Madricardo et al, 2017](https://doi.org/10.1038/sdata.2017.121) collected in the navigable channels of the Venice lagoon during a survey made in 2013 and reports depths relatives to the *Punta Salute level*. This dataset was coarsed by a factor 20 resulting in a coarse grid of resolution 10m x 10 m in order to reduce the processing time required by the interpolation procedures. In this file, coordinates are expressed in meters.

_see Bathymetry_scaling.pdf for details and visual explanation_
  
- **ADRI_CADEAU_NS** : this file introduces the gridded domain (resolution of 1/128 of degree) with coordinates expressed using longitude and latitude. Correspondent depth data have been neglected due to the coarser resolution with respect to the other files.


In order to produce a final bathymetry relative to the zero IGM 1942, the Corila2003 dataset is here used without depth correction while the depths of the coarsed2013 bathymetry are lowered by $24.3cm$ given that in 2013, the *Punta Saute level* was about $24.3cm$ lower than the *zero IGM 1942*.

## Required files
 - Bathymetry of the North Adriatic Sea (netCDF format) on which the Venice lagoon bathymetry is innested, in the following referred as _original grid_ (in this example _bathy_ADRI_CADEAU_NS.nc_)
 - Files (csv format) with depth data to update the original grid, in the following referred as _external files_ (in this example _Bathy_2003CORILA.csv_ and _Bathy_2013_coarsed.csv_)

## Required packages

- pandas
- geopandas
- xarray 
- matplotlib
- numpy 
- shapely.geometry
- math
- os
- typing

## _average_bathymetry.py_
The procedure to construct the Venice lagoon updated bathymetry consists of 5 steps: __Pre-processing__ data, __Averaging__ data recorded in the _external files_, __Weighing__ data given their occurrence on the grid-cells, __Cleaning__ of critical points (clusters of isolated points, land points, boundary points), __Adding__ the averaged-weighted bathymetry to the original bthymetry, __Exporting__ the results in different formats.

To start the procedure, the user has to modify the path where the necessary files are stored _fileLoc_ and the one where the produced files have to be saved _fileDest_: (they are both defined as the current folder by default).

In the following each step is described in details:

### PRE-PROCESSING 
The first step fulfilled is the conversion of the reference system of coordinates of the _external files_ in order to be comparable with the _original grid_ reference. To do so, all the _external files_ are converted into GeoDataFrame to exploit the built-in function that automatically convert the coordinate reference system (to_crs()). The initial and final reference system have to be specified changing the values of the variables _init_ref_ representing the initial reference system and _fin_ref_ that is the final one

>init_ref = "3004"
>
>fin_ref = "4326"

(_See the final remarks for more details about data_)

### AVERAGING
The first step of the algorithm used for the averaging consists of the construction of a more refined grid in which each cell of the _original grid_ (matrices defined on these macro-cells are identified in the source code as `bins` or `global`) is sub-divided into _nRef <sup>2</sup>_ sub-cells  (matrices defined on this refined grid are identified in the source code as `ref_`), where _nRef_ corresponds to the number of refinements desired and is a customizable parameter.

> nRef = 7 
> 
> This value depends on the density of samplings contained in the unstructured files, it is the largest value to certainly capture external data in each sub-cells (in area where data are present).

In each cell, two different histograms of dimensions $nRef \times nRef$ are constructed :
- the first histogram (`ref_occ`) records per each sub-cell (among the $nRef \times nRef$ sub-cells) the number of observations that fall within this sub-cell area, i.e. the number of occurences
- the second histogram (`ref_sum`) computes the sum of the depths of the observations that fall within this sub-cell area. 

At this point, a refined average depth is computed within each sub-cell dividing the value of the sum of the depths by the number of observation per each sub-cell. Sub-cells containing no-observations are considered land and set to $depth=0$. This could be done because $nRef$ was accurately chosen to insure the maximal possible refinement while keeping sub-cells dimensions larger than the resolution of the bathymetry observational datasets so that all sub-cells in water regions contain at least one observational data.

Once this computation is performed, for both datasets, the final sub-cells depth matrix is set equal to the depth in the sub-cells obtained with the first dataset (2003), then eventually erased by the depths of the sub-cells obtained with the second dataset (2016) where those sub-cells have depth values > 0. 
Finally, the average depth is computed within each macro-cell of the _original grid_ as a further average operation between all the sub-cells of this macro-cell (sum of the average depth in each sub-cell / $nRef ^2$) and the resulting value is assigned to the original macro-cell. The output of this operation is a matrix with the same dimension as the original grid containing in each cell the averaged value of the depths coming from the two combined average procedures.

### WEIGHING 
The weighing operation is performed to take into account the distribution of data in the decision of the caracterization of each cell (water or land), the depth values computed in the averaging step is weighed with respect to the percentage of occurrence of water into each cell. In order to do so, in each cell and for each _external file_ a matrix of boolean representing the occurreces of data considering every sub-cell is constructed (_True_ values means that at least one values falls into the cell, _False_ is assigned in the opposite case) and the union of these matrices with respect to the different _external file_ is computed. The resulting matrix is used to compute the percentage of water covering each cell. If the value computed is under a certain threshold, the cell is converted into land (depth equal to zero). In this way, all the cells that have only a small area covered by water are neglected into the reconstruction of the bathymetry (this allows to avoid enlargement of the surface covered by water through very shallow areas)

### CLEANING
The cleaning procedure consists in the visual correction of the outcome of the previous operations in order to have a result consistent with the reality.
Each cell is categorized into a categorical variable as follows:

| Value     |Present/Removed | Reason  |
| :---      |       :----:       |    :----    |
| **0**     | P       | Water cell |
| **1**     | R       | Not significant |
| **2**     | R       | Water coverage < 20% |
| **3**     | R       | Depth value given     |
| **4**     | P       | Depth value received   |
| **5**     | P       | Depth value manually changed      |
| **6**     | R       | Superposition with _original grid_ data of the North Adriatic Sea  |
| **7**     | R       | Depth value less than threshold      |
| **8**     | R       | Closure of the domain     |
| **9**     | P       | Depth value come from a re-distribution and averaged of adiacent cells     |


This operation is ad hoc for each case and is performed introducing 
>__operation__.loc[operation.latitude[_latitude index_],operation.longitude[_longitude index_]] = _correspondent value_

The same modifications has to be done on the __bathy__ variable

<ins>Remark:</ins> values 3, 4 and 9 have been introduced to guarantee the conservation of the volume of water. For values 3 and 4, this is done by adding the value of the depth of some particular cells to the one adiacent, based on the necessity to avoid the connection of cells or to represent land cells without neglecting important water channels. Regarding value 9 instead, the correspondent cells are deepened in order to conserve channel structures. This is done by averaging the amount of water contained in eliminated cells and redistributing it among cells with value 9.


### ADDING
Once the final bathymetry is ready, it is added to the _original grid_ of the North Adriatic Sea given at the beginning through the simple superposition of the values present in the two different files, the _original grid_ and the _lagoon bathymetry_

### EXPORTING
The produced files are exported in the .nc and .bin format


### Visualization
Open Test_case/Venice_avg_cleaned_bathymetry.qgz_ to visualize the results produced for the analysis performed.

