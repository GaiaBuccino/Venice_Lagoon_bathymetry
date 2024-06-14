# Venice Lagoon Bathymetry
Algorithm that updates the structured bathymetry of the North Adriatic Sea (resolution at 1/128 degree) adding the bathymetry of the Venice lagoon constructed by averaging values coming from more recent unstructured bathymetries (CORILA2003 and coarsed2013)

## Required files
 - Bathymetry of the North Adriatic Sea (netCDF format) on which the Venice lagoon bathymetry is innested, in the following referred as __original bathymetry__ (in the example _bathy_ADRI_CADEAU_NS.nc_)
 - Files (csv format) with depth data to update the original bathymetry, in the following referred as __external files__ (in the example _Bathy_2003CORILA.csv_ and _Bathy_2013_coarsed.csv_)

## _average_bathymetry.py_
The procedure to construct the Venice lagoon updated bathymetry consists of 5 steps: __Pre-processing__ data, __Averaging__ data recorded in the _external files_, __Weighing__ data given their occurrence on the grid-cells, __CLEANING__ OF THE ISOLATED SPOTS, __CONSERVATION__ OF THE VOLUME OF WATER, __Adding__ the averaged-weighted bathymetry to the original bthymetry, __Exporting__ the results in different formats.

To start the procedure, the user has to modify the path where the necessary files are stored _fileLoc_ and the one where the produced files have to be saved _fileDest_:
```
fileLoc = '...'
fileDest = '...' 
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

### Weighing 
the weighing operation is performed to take into account the distribution of data in the assignment of the depth value to an entire cell, the depth computed in the averaging step is weighed with respect to the percentage of occurrence of water into each cell. In order to do so, in each cell and for each _external file_ a matrix of boolean representing the occurreces of data considering every sub-cell is constructed (_True_ values means that at least one values falls into the cell, _False_ is assigned in the opposite case) and the union of these matrices with respect to the different _external file_ is computed. The resulting matrix is used to compute the percentage of water covering each cell. If the value computed is under a certain threshold, the cell is converted into land (depth equal to zero). In this way, all the cells that have only a small area covered by water are neglected into the reconstruction of the bathymetry (this allows to avoid enlargement of the surface covered by water through very shallow areas)

###

###

### Adding
Once the final bathymetry is ready, it is added to the _original bathymetry_ of the North Adriatic Sea given at the beginning through the simple superposition of the values present in the two different files, the _original bathymetry_ and the _lagoon bathymetry_

### Exporting
The produced files are exported in the .nc and .bin format

## Remark: External files


La profondità della batimetria della nuova griglia Venlag 64 si riferisce allo ZERO IGM (1942). Le profondità degli elementi triangolari della griglia sono stati calcolati attraverso la procedura seguente:
1.) il dataset Corila 2003 (che si riferisce allo ZERO IGM di 1942) è stato utilizzato per interpolare linearmente le profondità sui nodi di una griglia regolare di risoluzione 20 metri (creando un dataset di profondità che chiamiamo Corila-2003-gridded) che copre tutto il dominio di calcolo.
2.) ad ogni elemento triangolare della griglia Venlag 64 è stato assegnato la profondità media dei punti della griglia regolare (Corila-2003-gridded) contenuti nell'elemento triangolare della griglia Venlag 64.
3.) Il dataset 2013 (DOI: 10.1594/IEDA/323605) contiene dati di profondità di altissima risoluzione nei canali della laguna, si riferisce per le profondità al livello Punta Salute. Perciò, le profondità di questo dataset sono state aumentate di 23.5 cm per tenere in conto la differenza di profondità tra il Livello Punta Salute e lo zero IGM 1942.
4.) Poi, solo un punto ogni 5 punti in ogni direzione lat/lon del dataset 2013 (relativo allo zero IGM 1942) sono stati estratti (utilizzare la totalità dei dati avrebbe richiesto troppo tempo di calcolo) per creare una griglia intermedia che possiamo chiamare Coarsed-2013-Bati-IGM.

5.) In seguito, ogni elemento triangolare della griglia Venlag 64 per il quale almeno 75% della propria superficia era coperta dai dati della griglia Coarsed-2013-Bati-IGM, la profondità è stata aggiornata per essere uguale alla profondità media dei punti della griglia Coarsed-2013-Bati-IGM contenuti nell'elemento triangolare della griglia Venlag 64. 
6) Dove meno di 75% della superficie dell'elemento triangolare era sovrapposto con il dataset Coarsed-2013-Bati-IGM, una media ponderata tra i dati Corila-2003-gridded e i dati Coarsed-2013-Bati-IGM è stata calcolata.
7) Dove nessun dato di batimetria dei dataset 2003 o 2013 era disponibile, la profondità è stata interpolata linearmente dai dati di profondità della griglia Shelflag.

