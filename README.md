# bo_masks package
This package was developped in the need to construct masks for the DC2 catalog in order to mask the region where bright stars high luminosity might affect detected objects/images.
## It is subdivised in 3 classes :
### 1) call_dc2
Used to generate clean samples using some basic quality flags. One can use his own flags but "Default" ones are avaible depending on which DC2 catalog you are using.
### 2) radius_study
As we only want to mask regions affected by the luminosity of the star, we mask a circle around each star within a radius defined as the radius for which the relative density of clusters in the circle is inferior to a given value (modulable). Here we use 0.9.
### 3) generate_masks
Here we use the results from 2) and the catalog from 1) to generate a healsparse map which contains the masks. One can chose to save them in a healpix format which is not optimized for high resolution maps.

##Requirements:
Healsparse v1.9.0
