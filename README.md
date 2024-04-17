[![DOI](https://zenodo.org/badge/759888613.svg)](https://zenodo.org/doi/10.5281/zenodo.10986498)

# DC2 stellar masks package
This package was developped in the need to construct masks that covers regions where high luminosity of stars affects detected objects and the images in the DC2 context.
This work is highly based on Zilong Du's work during [his thesis](https://escholarship.org/uc/item/0tk4181k), but also on [HSC SSP article](https://academic.oup.com/pasj/article/doi/10.1093/pasj/psx047/4004646).

## It is subdivided in 3 modules :
### I) call_dc2
It opens the sample you want based on 
- the name of the catalog you ask for (must be a Catalog avaible in GCR except for the bin_star method which only require a catalog).
- some selection criterias of your will : general criteria applied to every catalog = conditions, criteria applied to 'galaxies' = conditions_galaxies, criteria applied to 'stars' = conditions_stars
### II) radius_study
As we only want to mask regions affected by the luminosity of the star, we mask a circle around each star within a radius defined as the radius for which the density of clusters in the circle over the mean density of the tract is inferior to a given value called critical_density. This value is customizable in the config file. 
### III) generate_masks
Here we use the results from 2) and the catalog from 1) to generate a healsparse map which contains the masks. One can chose to save them in a healpix format which is not optimized for high resolution maps.

## Note:
Future version of the module might allow one to perform masks on catalogs which are not in GCR

## Requirements:
Healsparse v1.9.0
GCRCatalogs v1.4.0


