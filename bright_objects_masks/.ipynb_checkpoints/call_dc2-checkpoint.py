import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter#, sample_filter
from GCRCatalogs import GCRQuery
from astropy.table import Table
from configparser import ConfigParser
import sys
"""
Errors gestion to be implemented
put API to describe funcs and variable roles
"""
class openDC2:
    
    def __init__(self, name='dc2_object_run2.2i_dr6_v2_with_addons_v2'):
        """__init__

        Parameters
        ----------
        name : str, optional
            name of the GCRCatalog to use, by default 'dc2_object_run2.2i_dr6_v2_with_addons_v2'

        Returns
        -------
        None
        """        
        self.name = name
        self.catalog = GCRCatalogs.load_catalog(name)
        return None

    def _flags(self, conditions=None):
        """_flags internal function used to treat the conditions to apply to GCRCatalog before Query

        Parameters
        ----------
        conditions : list of strings, optional
            if None : basic quality cuts are used, by default None

        Returns
        -------
        list of strings
            conditions being applied in _query
        """        
        if conditions == None:
            print("Default flags selected")
            filters = ["detect_isPrimary==True", "modelfit_CModel_flag_badCentroid==False", "base_SdssCentroid_flag==False", "base_PixelFlags_flag_edge==False"]
            filters += ["base_PixelFlags_flag_interpolatedCenter==False", "base_PixelFlags_flag_saturatedCenter==False", "base_PixelFlags_flag_bad==False"]
            filters += ["base_PixelFlags_flag_suspectCenter==False", "deblend_skipped==False", "base_PsfFlux_flag==False", "base_SdssShape_flag_psf==False"]
            filters += ["modelfit_DoubleShapeletPsfApprox_flag==False", "base_Blendedness_abs<=0.42169650342", "base_ClassificationExtendedness_flag==False", "snr_i_cModel > 5"]
            mag_bands = ["u", "g", "r", "i", "z", "y"]
            for mag in mag_bands:
                filters += [f"{mag}_base_PsfFlux_flag==False", f"{mag}_base_PixelFlags_flag_edge==False", f"{mag}_base_PixelFlags_flag_saturatedCenter==False"]
                if mag=="g" or mag=="r" or mag=="i":
                    filters += [f"{mag}_base_ClassificationExtendedness_flag==False"]
            mag_filters_ = mag_bands.copy().pop(3)
            for i in range(len(mag_filters_)): #We exclude i band from couples
                for j in range(i + 1, len(mag_filters_)):
                    filters += [f"(snr_{mag_filters_[i]}_cModel>3) | (snr_{mag_filters_[j]}_cModel>3)"]
            conditions = filters #Add all flags
        else :
            conditions = conditions
        return conditions
    
    def _query(self, conditions=None):
        """_query apply cuts on opened GCRCatalog

        Parameters
        ----------
        conditions : list of strings, optional
            if None : basic quality cuts are used, by default None

        Returns
        -------
        GCRQuery
            GCRQuery applied in catalog.get_quantities
        """        
        conditions_list = self._flags(conditions)
        if type(conditions_list) == list:
            if len(conditions_list)>1:
                filters = GCRQuery()
                for condition in conditions_list:
                    if type(condition) != str or len(condition) < 2: #if len cdt<2 there might be an error in the defined conditions (happened)
                        print("ERROR : condition must be of string type ex : detect_isPrimary==True") # do a verif str func ? 
                        break
                    filters &= GCRQuery(condition)
            else :
                filters = GCRQuery(conditions_list[0])
        elif type(conditions_list) == str:
            filters = GCRQuery(conditions_list)
        else : 
            print("ERROR condition must be of string type ex : detect_isPrimary==True or a list of condition in string type")
        return filters
    
    def open_cat(self, quantities = ['ra', 'dec', 'mag_i_cModel'], conditions = None, tract_list=None):
        """open_cat use get_quantities method for the selected GCRcatalog

        Parameters
        ----------
        quantities : list of strings, optional
            Catalog's parametters to get, by default ['ra', 'dec', 'mag_i_cModel']
        conditions : list of strings, optional
            General cuts to use (so called basic quality cut), by default None
        tract_list : list, optional
            Tract(s) to open the catalog in. If None : opens the full catalog, by default None

        Returns
        -------
        astropy Table
            Table of catalog objects with asked quantities
        """        
        if 'None' in str(conditions) and len(str(conditions).split(','))>1:
            filters = self._query()
            conditions.pop(conditions.index('None'))
            filters &= self._query(conditions = conditions)
        else :
            filters = self._query(conditions = conditions)
        if tract_list is not None:
            dc2 = self.catalog.get_quantities(quantities, native_filters=[tract_filter(tract_list)], filters=filters)
            print(f"DC2 catalog loaded (tract = {tract_list}, with filters)")
        else:
            dc2 = self.catalog.get_quantities(quantities, filters=filters)
            print("Full DC2 catalog loaded (with filters)")
        return Table(dc2) # Table format is better for future operations
    
    def galaxies(self, quantities=['ra','dec','mag_i_cModel'], conditions=None, conditions1 = ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"], tract_list=None):
        """get galaxies from the catalog

        Parameters
        ----------
        quantities : list, optional
            Catalog's parametters to get, by default ['ra','dec','mag_i_cModel']
        conditions : list, optional
            General cuts to use (so called basic quality cut), by default None = basic quality cuts
        conditions1 : list, optional
            Cuts to use to select galaxies, by default ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"]
        tract_list : int or list, optional
            Tract(s) to open the catalog in. by default None : opens the full catalog

        Returns
        -------
        astropy Table
            Catalog of galaxies
        """
        if conditions is None:
            conditions = [str(conditions)]
        dc2_galaxies = self.open_cat(quantities, list(conditions) + list(conditions1), tract_list)
        return dc2_galaxies
    
    def stars(self, quantities=['ra','dec','mag_i_cModel'], conditions=None, conditions1=["extendedness==0"], tract_list=None):
        """Isolate stars from the 'cleaned' catalog

        Parameters
        ----------
        quantities : list, optional
            Catalog's parametters to get, by default ['ra','dec','mag_i_cModel']
        conditions : list, optional
            General cuts to use (so called basic quality cut), by default None
        conditions1 : list, optional
            Cuts to use to select stars, by default ["extendedness==0"]
        tract_list : int or list, optional
            Tract(s) to open the catalog in. If None : opens the full catalog, by default None

        Returns
        -------
        astropy Table
            Catalog of stars
        """      
        if conditions is None:
            conditions = [str(conditions)]
        dc2_stars = self.open_cat(quantities, list(conditions) + list(conditions1), tract_list)
        return dc2_stars
    
    def bin_cat(self, catalog, quantities="mag_i_cModel", bins=[0, 17, 18, 20, 22, 24]):
        """bin_cat sorts objects from a given catalog in bins of a selected quantity. Often used to bin stars by magnitude

        Parameters
        ----------
        catalog : astropy Table
            Any type of objects from the initial GCRCatalog.
        quantities : str, list of str, optional
            Quantities to bin the catalog in, by default "mag_i_cModel"
        bins : list, optional
            Bins of 'quantities', by default [0, 17, 18, 20, 22, 24]

        Returns
        -------
        list
            List containing each catalog's bins
        """        
        binned_cat = []
        for i in range(len(bins)-1):
            v_min = bins[i]
            v_max = bins[i+1]
            binned_cat.append(catalog[(catalog[quantities]>v_min) & (catalog[quantities]<v_max)])
        return binned_cat
    
    def galaxies_with_neighbours_tracts(self, quantities=['ra','dec','mag_i_cModel'], conditions=None, conditions1=["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"], tract_list = None):
        """galaxies_with_neighbours_tracts opens all neighbours tracts catalog (often galaxies) for a given tract or several tracts

        Parameters
        ----------
        quantities : list, optional
            Catalog's parametters to get, by default ['ra','dec','mag_i_cModel']
        conditions : str, list of str, optional
            General cuts to use (so called basic quality cut), by default None
        conditions1 : str, list of str, optional
            Cuts to use to select galaxies, by default ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"]
        tract_list : int, list, optional
            Tract(s) to open the catalog in, by default None

        Returns
        -------
        astropy Table
            Astropy Table containing selected tract + neighbours galaxies
        """        
        neighbour_tracts = Table.read("/sps/lsst/groups/clusters/amico_validation_project/catalogs/DC2/dc2_neighbours.fits")
        neighbour_list = neighbour_tracts["list_of_neighbour_tiles"][neighbour_tracts["tile"]==int(tract_list)][0].split(',')
        dc2_galaxies = self.galaxies(quantities, conditions, conditions1, tract_list=neighbour_list)
        return dc2_galaxies