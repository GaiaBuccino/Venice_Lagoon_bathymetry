# Venice Lagoon Bathymetry
Algorithm that updates the structured bathymetry of the North Adriatic Sea (resolution at 1/128 degree) adding the bathymetry of the Venice lagoon constructed by averaging values coming from more recent unstructured bathymetries (CORILA2003 and coarsed2013)

## Required files
 - Bathymetry of the North Adriatic Sea (netCDF format) on which the Venice lagoon bathymetry is innested, in the following referred as __original bathymetry__ (in the example _bathy_ADRI_CADEAU_NS.nc_)
 - Files (csv format) with depth data to update the original bathymetry, in the following referred as __external files__ ( in the example _Bathy_2003CORILA.csv_ and _Bathy_2013_coarsed.csv_)
