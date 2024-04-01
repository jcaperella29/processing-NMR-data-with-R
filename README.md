# processing-NMR-data-with-R
A script for processing  NMR data using R. The script outputs a csv with the entire processed spectra as well as a text file containing some of peeks .The later can be used to map the peeks to chemicals.
First needed libraries are imported and the samples are extracted.
then the data is processed  using the  "Rnmr1D" package 
the spectra is then ploted using a stacked plot
Then normlization and clustering is performed. The clustering is by concentration in ppm.
The clustering results are  then plotted.
Next PCA is performed 
Finally resutls are exported
