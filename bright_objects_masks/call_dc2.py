import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter  # , sample_filter
from GCRCatalogs import GCRQuery
from astropy.table import Table, vstack
from configparser import ConfigParser
import sys
from lsst.daf.butler import Butler

"""
Errors gestion to be implemented.
"""


class OpenDC2:

    def __init__(self, name="dc2_object_run2.2i_dr6_v2_with_addons_v2", butler = None, collection = None):
        """__init__

        Parameters
        ----------
        name : str, optional
            name of the GCRCatalog to use, by default 'dc2_object_run2.2i_dr6_v2_with_addons_v2'

        Returns
        -------
        None
        """
        if butler is None:
            self.name = name
            self.catalog = GCRCatalogs.load_catalog(self.name)
            self.butler = butler
        else :
            self.butler = Butler(butler, collections=collection)
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
        if conditions is None:
            print("Default flags selected")
            filters = [
                "detect_isPrimary==True",
                "modelfit_CModel_flag_badCentroid==False",
                "base_SdssCentroid_flag==False",
                "base_PixelFlags_flag_edge==False",
            ]
            filters += [
                "base_PixelFlags_flag_interpolatedCenter==False",
                "base_PixelFlags_flag_saturatedCenter==False",
                "base_PixelFlags_flag_bad==False",
            ]
            filters += [
                "base_PixelFlags_flag_suspectCenter==False",
                "deblend_skipped==False",
                "base_PsfFlux_flag==False",
                "base_SdssShape_flag_psf==False",
            ]
            filters += [
                "modelfit_DoubleShapeletPsfApprox_flag==False",
                "base_Blendedness_abs<=0.42169650342",
                "base_ClassificationExtendedness_flag==False",
                "snr_i_cModel > 5",
            ]
            mag_bands = ["u", "g", "r", "i", "z", "y"]
            for mag in mag_bands:
                filters += [
                    f"{mag}_base_PsfFlux_flag==False",
                    f"{mag}_base_PixelFlags_flag_edge==False",
                    f"{mag}_base_PixelFlags_flag_saturatedCenter==False",
                ]
                if mag == "g" or mag == "r" or mag == "i":
                    filters += [f"{mag}_base_ClassificationExtendedness_flag==False"]
            mag_filters_ = mag_bands.copy().pop(3)
            for i in range(len(mag_filters_)):  # We exclude i band from couples
                for j in range(i + 1, len(mag_filters_)):
                    filters += [
                        f"(snr_{mag_filters_[i]}_cModel>3) | (snr_{mag_filters_[j]}_cModel>3)"
                    ]
            conditions = filters  # Add all flags
        else:
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
            if len(conditions_list) > 1:
                filters = GCRQuery()
                for condition in conditions_list:
                    if (
                        type(condition) != str or len(condition) < 2
                    ):  # if len cdt<2 there might be an error in the defined conditions (happened)
                        print(
                            "ERROR : condition must be of string type ex : detect_isPrimary==True"
                        )  # do a verif str func ?
                        break
                    filters &= GCRQuery(condition)
            else:
                filters = GCRQuery(conditions_list[0])
        elif type(conditions_list) == str:
            filters = GCRQuery(conditions_list)
        else:
            print(
                "ERROR condition must be of string type ex : detect_isPrimary==True or a list of condition in string type"
            )
        return filters

    def open_cat(
        self, quantities=["ra", "dec", "mag_i_cModel"], conditions=None, tract_list=None, datasetType=None
    ):
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
        if self.butler is None:
            if "None" in str(conditions) and len(str(conditions).split(",")) > 1:
                filters = self._query()
                conditions.pop(conditions.index("None"))
                filters &= self._query(conditions=conditions)
            else:
                filters = self._query(conditions=conditions)
            if tract_list is not None:
                dc2 = self.catalog.get_quantities(
                    quantities, native_filters=[tract_filter(tract_list)], filters=filters
                )
                print(f"DC2 catalog loaded (tract = {tract_list}, with filters)")
            else:
                dc2 = self.catalog.get_quantities(quantities, filters=filters)
                print("Full DC2 catalog loaded (with filters)")
            return Table(dc2)
        else :
            if tract_list is not None:
                registry = self.butler.registry
                tract_cdt = ""
                for tract in tract_list:
                    tract_cdt += str(tract) + ","
                tract_cdt = "(" + tract_cdt[:-1] + ")"
                tract_query = " and tract in " + tract_cdt
                datasets = list(registry.queryDatasets(datasetType, where=conditions + tract_query))
                if datasetType == "the_monster_20250219":
                    if len(datasets) > 1:
                        dc2 = self.butler.get(datasets[0], parameters={"columns": quantities}).asAstropy()
                        for i in range(len(datasets) - 1):
                            dc2 = vstack([dc2, self.butler.get(datasets[i + 1], parameters={"columns": quantities}).asAstropy()])
                    else:
                        dc2 = self.butler.get(datasets[0], parameters={"columns": quantities}).asAstropy()
                else :
                    if len(datasets) > 1:
                        dc2 = self.butler.get(datasets[0], parameters={"columns": quantities})
                        for i in range(len(datasets) - 1):
                            dc2 = vstack([dc2, self.butler.get(datasets[i + 1], parameters={"columns": quantities})])
                    else:
                        dc2 = self.butler.get(datasets[0], parameters={"columns": quantities})
            else:
                print("Error : tract_list must be specified when using butler")
            return dc2

    def galaxies(
        self,
        quantities=["ra", "dec", "mag_i_cModel"],
        conditions=None,
        conditions1=["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"],
        tract_list=None,
        datasetType=None
    ):
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
        if self.butler is not None:
            if conditions is None:
                print("No conditions specified, need to use at least one in butler configuration.")
            else :
                dc2_galaxies = self.open_cat(
                    quantities, conditions, tract_list, datasetType=datasetType
                )
            if len(conditions1)==0:
                print("No conditions1 specified, returning catalog with conditions only.")
                return dc2_galaxies
            else: #put a warning later around here
                for condition in conditions1:
                    parts = condition.split("==")
                    if len(parts) == 2:
                        column_name = parts[0]
                        value = parts[1]
                        if column_name in dc2_galaxies.colnames:
                            dc2_galaxies = dc2_galaxies[dc2_galaxies[column_name] == value]

                    elif ">" in condition:
                        parts = condition.split(">")
                        column_name = parts[0]
                        value = float(parts[1])
                        if column_name in dc2_galaxies.colnames:
                            dc2_galaxies = dc2_galaxies[dc2_galaxies[column_name] > value]
                
                    elif "<" in condition:
                        parts = condition.split("<")
                        column_name = parts[0]
                        value = float(parts[1])
                        if column_name in dc2_galaxies.colnames:
                            dc2_galaxies = dc2_galaxies[dc2_galaxies[column_name] < value]
        else :
            if conditions is None:
                conditions = [str(conditions)]
            dc2_galaxies = self.open_cat(
                quantities, list(conditions) + list(conditions1), tract_list
            )
        return dc2_galaxies

    def stars(
        self,
        quantities=["ra", "dec", "mag_i_cModel"],
        conditions=None,
        conditions1=["extendedness==0"],
        tract_list=None,
        datasetType=None
    ):
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
        if self.butler is not None:
            if conditions is None:
                dc2_stars = self.open_cat(
                    quantities, conditions1, tract_list, datasetType=datasetType
                )
            if len(conditions1)==0:
                print("No conditions1 specified, returning catalog with conditions only.")
                return dc2_stars
            else: #put a warning later around here
                for condition in conditions1:
                    parts = condition.split("==")
                    if len(parts) == 2:
                        column_name = parts[0]
                        value = parts[1]
                        if column_name in dc2_stars.colnames:
                            dc2_stars = dc2_stars[dc2_stars[column_name] == value]

                    elif ">" in condition:
                        parts = condition.split(">")
                        column_name = parts[0]
                        value = float(parts[1])
                        if column_name in dc2_stars.colnames:
                            dc2_stars = dc2_stars[dc2_stars[column_name] > value]

                    elif "<" in condition:
                        parts = condition.split("<")
                        column_name = parts[0]
                        value = float(parts[1])
                        if column_name in dc2_stars.colnames:
                            dc2_stars = dc2_stars[dc2_stars[column_name] < value]
        else :
            if conditions is None:
                conditions = [str(conditions)]
            dc2_stars = self.open_cat(
                quantities, list(conditions) + list(conditions1), tract_list
            )
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
        try:
            for i in range(len(bins) - 1):

                v_min = float(bins[i])
                v_max = float(bins[i + 1])
                binned_cat.append(
                    catalog[
                        (catalog[quantities] > v_min) & (catalog[quantities] < v_max)
                    ]
                )
        except KeyError:
            print(
                "Error: The specified column name 'quantities' does not exist in the catalog."
            )
        except ValueError:
            print("Error: Unable to convert bin values to float.")
        return binned_cat

    def galaxies_with_neighbours_tracts(
        self,
        quantities=["ra", "dec", "mag_i_cModel"],
        conditions=None,
        conditions1=["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"],
        tract_list=None,
        neighbours_tracts_path = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/DC2/dc2_neighbours.fits",
        datasetType = None
    ):
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
        neighbour_tracts = Table.read(
            neighbours_tracts_path
        )
        if len(tract_list) > 1 and type(tract_list) == list:
            neighbour_list = []
            neighbour_list.append(
                neighbour_tracts["list_of_neighbour_tiles"][
                    neighbour_tracts["tile"] == int(tract_list[0])
                ][0].split(",")
            )
            for i in range(len(tract_list) - 1):
                neighbour_list.append(
                    neighbour_tracts["list_of_neighbour_tiles"][
                        neighbour_tracts["tile"] == int(tract_list[i + 1])
                    ][0].split(",")
                )
                neighbour_list[i + 1] = neighbour_list[i] + list(
                    set(neighbour_list[i + 1]) - set(neighbour_list[i])
                )
            neighbour_list = neighbour_list[-1]
        elif type(tract_list) == list:
            neighbour_list = neighbour_tracts["list_of_neighbour_tiles"][
                neighbour_tracts["tile"] == int(tract_list[0])
            ][0].split(",")
        elif type(tract_list) == str or type(tract_list) == int:
            neighbour_list = neighbour_tracts["list_of_neighbour_tiles"][
                neighbour_tracts["tile"] == int(tract_list)
            ][0].split(",")
        dc2_galaxies = self.galaxies(
            quantities, conditions, conditions1, tract_list=neighbour_list, datasetType=datasetType
        )

        return dc2_galaxies
