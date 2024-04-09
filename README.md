# Venice_Lagoon_hydrodynamics

Repository to collect the comparison between hydrodynamics in the Adriatic Sea and the same area including the Venice Lagoon.

### Required file: 
>- Batimetria.dat


## 1. create_bathy_ADRI.m: produrre e scaricare i file

    - Bathymetry 792 x 424 Adriatic Sea + Lagoon (lat, lon, depth) [completa_1_128.bin]
    - Bathymetry 792 x 424 Lagoon + Adriatic Sea set to zero (lat, lon, depth) [lagoon_ADRIsize.bin]
    - Bathymetry 792 x 424 Adriatic Sea + Lagoon set to zero [ADRI_1_128_NO_Lagoon.bin]

## 2. file read_bathy.py: 

uplloads the .bin files coming from matlab and convert them into .csv
carica i file di batimetria unstructured e li trasforma nel sistema di riferimento EPSG:4326 (Coarsed2013) 

% NB: la profondità considerata è la prima registrata nei file di batimetrie

## 3. interpolation.py:

carica i file salvati da read_bathy.py
allarga la griglia di un'ampiezza pari a metà del passo di discretizzazione spaziale, i punti della griglia strutturata diventano così i centro-cella
interpola i valori delle profondità dei dati non strutturati sulla griglia strutturata (media dei punti unstructured in ogni cella structured) con hist2d
sostituisce ai valori di profondità del file Laguna + Adriatico a zero i valori mediati
si ottiene il file adjusted_bathy.csv e il plot Adjusted_bathymetry.png

## 4. bathy_bin_construction.py:

costruzione della batimetria finale come (batimetria Adriatico + Laguna a zero) + (batimetria adjusted)
