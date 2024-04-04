# Venice_Lagoon_hydrodynamics

Repository to collect the comparison between hydrodynamics in the Adriatic Sea and the same area including the Venice Lagoon.

1. ##DA MATLAB: produrre e scaricare i file

    - Batimetria 792 x 424 Adriatico + Laguna (lat, lon, depth) [completa_1_128.bin]
    - Batimetria 792 x 424 Laguna + Adriatico a zero (lat, lon, depth) [lagoon_ADRIsize.bin]
    - Batimetria 792 x 424 Adriatico + Laguna a zero [ADRI_1_128_NO_Lagoon.bin]

2. ##file read_bathy.py: 

carica i file .bin provenienti da matlab e li trasforma in csv
carica i file di batimetria unstructured e li trasforma nel sistema di riferimento EPSG:4326 (Coarsed2013) 

% NB: la profondità considerata è la prima registrata nei file di batimetrie

3. ##interpolation.py:

carica i file salvati da read_bathy.py
allarga la griglia di un'ampiezza pari a metà del passo di discretizzazione spaziale, i punti della griglia strutturata diventano così i centro-cella
interpola i valori delle profondità dei dati non strutturati sulla griglia strutturata (media dei punti unstructured in ogni cella structured) con hist2d
sostituisce ai valori di profondità del file Laguna + Adriatico a zero i valori mediati
si ottiene il file adjusted_bathy.csv e il plot Adjusted_bathymetry.png

4. ##bathy_bin_construction.py:

costruzione della batimetria finale come (batimetria Adriatico + Laguna a zero) + (batimetria adjusted)
