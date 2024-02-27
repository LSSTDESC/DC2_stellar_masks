import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import hpgeom as hpg
import healsparse as hsp
import healpy as hp
from astropy.io import fits
import radius_study
import call_dc2

class masks:
    def __init__(self, name='dc2_object_run2.2i_dr6_v2_with_addons_v2', theta_bins = np.logspace(np.log10(0.5),np.log10(50), 50), tract_list = None, 
                        quantities = ['ra','dec','mag_i_cModel'], conditions=None, conditions_galaxies=["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"], 
                        conditions_stars=["extendedness==0"], binned_quantity="mag_i_cModel", bins=[0, 17, 18, 20, 22, 24], density_ratio = None, critical_density = 0.9,
                        nside_coverage = 32, nside_sparse = 131072):
        """__init__

        Parameters
        ----------
        name : str, optional
            name of the GCRCatalog to use, by default 'dc2_object_run2.2i_dr6_v2_with_addons_v2'
        theta_bins : _type_, optional
            Radius in which you'll look for the density of galaxies around the star, by default np.logspace(np.log10(0.5),np.log10(50), 50)
        tract_list : _type_, optional
            Tract(s) to open the catalog in, by default None
        quantities : _type_, optional
            Catalog's parametters to get, by default  ['ra','dec','mag_i_cModel']
        conditions : _type_, optional
            General cuts to use (so called basic quality cut), by default None
        conditions_galaxies : _type_, optional
            Cuts to use to select galaxies, by default ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"]
        conditions_stars : _type_, optional
            Cuts to use to select stars, by default ["extendedness==0"]
        binned_quantity : str, optional
            Quantities to bin the catalog,, by default "mag_i_cModel"
        bins : list, optional
            Bins of 'quantities', by default [0, 17, 18, 20, 22, 24]
        density_ratio : _type_, optional
            density_ratio(r) for each bin of stars, by default None
        critical_density : float, optional
            critical density to get around each star, by default 0.9. Will have an impact on the cut radius around stars
        nside_coverage : int, optional
            Resolution parameter of coverage healsparse map, by default 32
        nside_sparse : int, optional
            Resolution parameter of sparse healsparse map, by default 131072
        """        
        self.name = name
        self.theta_bins = theta_bins
        self.tract_list = tract_list
        self.quantities = quantities
        self.conditions = conditions
        self.conditions_galaxies = conditions_galaxies
        self.conditions_stars = conditions_stars
        self.binned_quantity = binned_quantity
        self.bins = bins
        self.density_ratio = density_ratio
        self.critical_density = critical_density
        critical_radius = radius_study.critical_radius(name=self.name, theta_bins=self.theta_bins, tract_list=self.tract_list, quantities=self.quantities, conditions=self.conditions, conditions_galaxies=self.conditions_galaxies, conditions_stars=self.conditions_stars, binned_quantity=self.binned_quantity,bins=self.bins)
        self.radius = critical_radius.get_critical_radius(density_ratio = self.density_ratio, critical_density=self.critical_density)
        self.nside_coverage = nside_coverage
        self.nside_sparse = nside_sparse
        self.openDC2 = call_dc2.openDC2(name=self.name)
        if self.tract_list is None:
            self.galaxies = self.openDC2.galaxies(quantities=self.quantities, conditions=self.conditions, conditions1=self.conditions_galaxies)
            self.stars = self.openDC2.stars(quantities=self.quantities, conditions=self.conditions, conditions1=self.conditions_stars)
        else :
            self.galaxies = self.openDC2.galaxies(quantities=self.quantities, conditions=self.conditions, conditions1=self.conditions_galaxies, tract_list = tract_list)
            self.stars_cat = self.openDC2.stars(quantities=self.quantities, conditions=self.conditions, conditions1=self.conditions_stars, tract_list = tract_list)
        if binned_quantity is None:
            self.stars = [self.stars_cat]
        else : 
            self.stars = self.openDC2.bin_cat(self.stars_cat, quantities=self.binned_quantity, bins=self.bins)

    def create_healsparse_masks(self):
        """create_healsparse_masks Creates the healsparse map containing masks with previously given parameters

        Returns
        -------
        healsparse map
            Masks in a healsparse map format
        """        
        ra_min, ra_max = min(self.galaxies ["ra"]), max(self.galaxies ["ra"])
        dec_min, dec_max = min(self.galaxies ["dec"]), max(self.galaxies ["dec"])
        healsparse_masks = hsp.HealSparseMap.make_empty(self.nside_coverage, self.nside_sparse, np.bool_)
        polygon = hsp.Polygon(ra=[ra_min, ra_max, ra_max, ra_min],
                      dec=[dec_min, dec_min, dec_max, dec_max],
                      value=False)
        healsparse_masks &= polygon
        for i in range(len(self.stars)):
            stars = self.stars[i]
            radius = self.radius[i]
            print(len(stars))
            for star in stars :
                circle = hsp.Circle(ra=star["ra"], dec=star["dec"], radius=radius/3600, value=True) #radius value in degrees
                healsparse_masks |= circle
        return healsparse_masks

    def write_heaslparse_mask(self, healsparse_masks, outpath = None):
        """write_heaslparse_mask save created masks

        Parameters
        ----------
        healsparse_masks : healsparse map
            healsparse_map
        outpath : str, optional
            path in which to save masks, by default ''

        Returns
        -------
        """    
        if outpath is None:
            outpath = f'bo_masks_{self.tract_list}.hs'
        healsparse_masks.write(outpath , clobber=True)
        return None

    def convert_maks(self, healsparse_masks, nside_hp = int(4096*2), reduction_method = "sum"):
        """convert_maks uses healsparse module to convert healsparse map to healpix map

        Parameters
        ----------
        healsparse_masks : healsparse map
            healsparse map
        nside_hp : int, optional
            resolution of the output converted map, by default int(4096*2)
        reduction_method : str, optional
            reduction method from healsparse module, by default "sum"

        Returns
        -------
        ra, dec, map
            center positions of each healpix pixel and value of each pixel
        """        
        healpix_map = healsparse_masks.generate_healpix_map(nside= nside_hp, reduction = reduction_method)
        ra_hp, dec_hp = hp.pix2ang(nside_hp, np.arange(len(healpix_map)), lonlat = True, nest = True)
        ra_hp_s = ra_hp[(ra_hp<ra_max) & (ra_hp>ra_min) & (dec_hp<dec_max) & (dec_hp>dec_min)]
        dec_hp_s = dec_hp[(ra_hp<ra_max) & (ra_hp>ra_min) & (dec_hp<dec_max) & (dec_hp>dec_min)]
        hp_map_s = healpix_map[(ra_hp<ra_max) & (ra_hp>ra_min) & (dec_hp<dec_max) & (dec_hp>dec_min)]
        ra_hp, dec_hp, hp_map = ra_hp_s, dec_hp_s, hp_map_s
        return ra_hp, dec_hp, hp_map

    def write_healpix_masks(self, ra_hp, dec_hp, hp_map, outpath = None):
        """write_healpix_masks saves converted healpix map

        Parameters
        ----------
        ra_hp : array
            right ascension of each healpix pixel
        dec_hp : array
            declination of each healpix pixel
        hp_map : _type_
            value of each healpix pixel
        """        
        if outpath is None :
            outpath = f"bo_masks_{self.tract_list}.fits"
        hp_map = Table({"mask" : hp_map})
        hp_map["mask_bool"] = 1
        hp_map["mask_bool"][hp_map["mask"]/256 > 0.13] = 0 #To improve but says that >13% of surface covered = masked pixel
        print("Now saving...")
        hp_map["ra"], hp_map["dec"] = ra_hp, dec_hp
        header = fits.Header()
        header['RANGE_RA'] = f'{[min(ra_hp), max(ra_hp)]}'
        header['RANGE_DEC'] = f'{[min(dec_hp), max(dec_hp)]}'
        header['NSIDE'] = f'{nside_hp}'
        header['NEST'] = 'nest = True'
        header['DOWNGRADE_METHOD'] = 'sum for healsparse + > 13% of pixel masked'
        hdu = fits.BinTableHDU(hp_map, header=header)
        hdu.writeto(outpath, overwrite=True)
        return None