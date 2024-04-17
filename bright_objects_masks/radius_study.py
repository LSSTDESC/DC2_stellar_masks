import numpy as np
import sys, os
import call_dc2
import yaml

sys.path.append("_multiprocessing")
from _mprocessing.s_multiprocessing import Multiprocessing


class Calculate_density:
    """Class internally used for calculations"""

    def __init__(self, galaxies_tract, galaxies_tract_neighbours, stars):
        """__init__

        Parameters
        ----------
        galaxies_tract : astropy Table
            Galaxies of the tract in which computations will be done
        galaxies_tract_neighbours : astropy Table
            Galaxies from neighbouring tracts of the selected one
        stars : astropy Table
            Stars of the tract in which computations will be done, best practice is to use binned galaxies

        Returns
        -------
        """
        self.galaxies_tract = galaxies_tract
        self.galaxies_tract_neighbours = galaxies_tract_neighbours
        self.stars = stars
        return None

    ### Geometry
    def _distance(self, ra1, dec1, ra2, dec2):
        """_distance computes distance between two objects

        Parameters
        ----------
        ra1 : float or array
            right ascension of object 1 [deg]

        dec1 : float or array
            declination of object 1 [deg]

        ra2 : float or array
            right ascension of object 2 [deg]

        dec2 : float or array
            declination of object 2 [deg]

        Returns
        -------
        float or array
            distance between the two objects [arcseconds]
        """
        delta_ra = ra2 - ra1
        delta_dec = dec2 - dec1
        d1 = np.sqrt(delta_ra**2 + delta_dec**2)
        d1 = (d1) * 3600
        return d1

    def _surface_circle(self, theta, dec):
        """_surface_circle calculates the surface of a circle for a given theta radius with the center at dec

        Parameters
        ----------
        theta : float, array
            radius of the circle [arcsec]
        dec : float, array
            position of the center of the circle [deg]

        Returns
        -------
        float, array
            surface of the circle [arcsec**2]
        """
        return np.pi * (theta**2) * np.sin(np.radians(90 - dec))

    def _surface_trapeze(self):
        """_surface_trapeze computes the surface of one tract with the assumption that it's a trapeze

        Returns
        -------
        float
            area of the tract [arcsec**2]
        """
        tract_galaxies = self.galaxies_tract
        ra_galaxies, dec_galaxies = np.array(tract_galaxies["ra"]), np.array(
            tract_galaxies["dec"]
        )
        dec_min, dec_max = min(dec_galaxies), max(dec_galaxies)
        condition_B = (abs(dec_galaxies) >= abs(dec_min) - 0.2) & (
            abs(dec_galaxies) <= abs(dec_min) + 0.2
        )  # Find the small side
        condition_b = (abs(dec_galaxies) >= abs(dec_max) - 0.1) & (
            abs(dec_galaxies) <= abs(dec_max) + 0.1
        )  # Find the larger side
        ra_min_b, ra_max_b = min(ra_galaxies[condition_b]), max(
            ra_galaxies[condition_b]
        )
        ra_min_B, ra_max_B = min(ra_galaxies[condition_B]), max(
            ra_galaxies[condition_B]
        )
        small_side = ra_max_b - ra_min_b
        large_side = ra_max_B - ra_min_B
        mean_side = (large_side + small_side) / 2
        delta_dec = dec_max - dec_min
        center_dec = (
            dec_min + dec_max
        ) / 2  # Because we are on a sphere, surface needs to be corrected with dec position
        tract_surface = (
            mean_side * delta_dec * np.sin(np.radians(90 - center_dec)) * (3600**2)
        )
        return tract_surface

    def _tract_density(self):
        """_tract_density calculates the density of clusters in a tract

        Returns
        -------
        float
            density of the tract [arcsec**(-2)]
        """
        tract_galaxies = self.galaxies_tract
        galaxies_number = len(tract_galaxies["ra"])
        tract_density = galaxies_number / self._surface_trapeze()
        return tract_density

    def get_density_around_unique_star(
        self, theta_bins=np.logspace(np.log10(0.5), np.log10(50), 50)
    ):
        density = np.zeros((len(self.stars), len(theta_bins)))
        print("in get_density_around_unique_star", self.stars)
        print(self.stars["ra"])
        print(len(self.stars["ra"]))
        condition = np.array(
            [
                (
                    abs(self.galaxies_tract_neighbours["ra"] - self.stars["ra"][i])
                    * 3600
                    <= max(theta_bins)
                )
                | (
                    abs(self.galaxies_tract_neighbours["dec"] - self.stars["dec"][i])
                    * 3600
                    <= max(theta_bins)
                )
                for i in range(len(self.stars["ra"]))
            ]
        )
        for i, star in enumerate(self.stars):
            ra_galaxies, dec_galaxies = (
                self.galaxies_tract_neighbours["ra"][condition[i]],
                self.galaxies_tract_neighbours["dec"][condition[i]],
            )  # Restrain number of galaxies for which we compute distance
            distance = self._distance(
                ra_galaxies, dec_galaxies, star["ra"], star["dec"]
            )  # Distance between star and galaxies within ra or dec max theta
            number_galaxies = np.array(
                [
                    len(distance[distance <= theta_bins[j]])
                    for j in range(len(theta_bins))
                ]
            )  # Number of objects within distance theta
            surface = self._surface_circle(theta=theta_bins, dec=star["dec"])
            density[i] = (
                number_galaxies / surface
            )  # true density for each theta_bin will then be density / len(stars)
        return density

    def get_density_around_stars(
        self, theta_bins=np.logspace(np.log10(0.5), np.log10(50), 50)
    ):
        """get_density_around_stars gives density(r) around bright stars in one tract, good practice is to give it magnitude binned stars

        Parameters
        ----------
        theta_bins : array, optional
            Radius in which you'll look for the density of galaxies around the star, by default np.logspace(np.log10(0.5),np.log10(50), 50)

        Returns
        -------
        array
            density as a function of theta for each bin of stars if stars are binned
        """
        density = np.zeros(len(theta_bins))
        condition = np.array(
            [
                (
                    abs(self.galaxies_tract_neighbours["ra"] - self.stars["ra"][i])
                    * 3600
                    <= max(theta_bins)
                )
                | (
                    abs(self.galaxies_tract_neighbours["dec"] - self.stars["dec"][i])
                    * 3600
                    <= max(theta_bins)
                )
                for i in range(len(self.stars["ra"]))
            ]
        )
        sum_density = 0
        for i, star in enumerate(self.stars):
            ra_galaxies, dec_galaxies = (
                self.galaxies_tract_neighbours["ra"][condition[i]],
                self.galaxies_tract_neighbours["dec"][condition[i]],
            )  # Restrain number of galaxies for which we compute distance
            distance = self._distance(
                ra_galaxies, dec_galaxies, star["ra"], star["dec"]
            )  # Distance between star and galaxies within ra or dec max theta
            number_galaxies = np.array(
                [
                    len(distance[distance <= theta_bins[j]])
                    for j in range(len(theta_bins))
                ]
            )  # Number of objects within distance theta
            surface = self._surface_circle(theta=theta_bins, dec=star["dec"])
            sum_density += (
                number_galaxies / surface
            )  # true density for each theta_bin will then be density / len(stars)
        number_stars = len(self.stars["ra"])
        if number_stars != 0:
            density = sum_density / number_stars  # Mean density over stars
        else:
            print("Error : division by zero encountered, setting density to nan")
            density = np.nan
        return density


class Critical_radius:
    def __init__(
        self, config_file=os.path.dirname(__file__) + "/../config/config.yaml"
    ):
        """self,
            name="dc2_object_run2.2i_dr6_v2_with_addons_v2",
            theta_bins=np.logspace(np.log10(0.5), np.log10(50), 50),
            tract_list=None,
            quantities=["ra", "dec", "mag_i_cModel"],
            conditions=None,
            conditions_galaxies=["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"],
            conditions_stars=["extendedness==0"],
            binned_quantity="mag_i_cModel",
            bins=[0, 17, 18, 20, 22, 24],
        ):"""
        """__init__ Here we look for the radius of the circle area we'll mask around each bright star

        Parameters
        ----------
        name : str, optional
            name of the GCRCatalog to use, by default 'dc2_object_run2.2i_dr6_v2_with_addons_v2'
        theta_bins : array, optional
            Radius in which you'll look for the density of galaxies around the star, by default np.logspace(np.log10(0.5),np.log10(50), 50)
        tract_list : list, optional
            Tract(s) to open the catalog in. by default None
        quantities : list, optional
            Catalog's parametters to get, by default ['ra','dec','mag_i_cModel']
        conditions : list, optional
            General cuts to use (so called basic quality cut), by default None
        conditions_galaxies : list, optional
            Cuts to use to select galaxies, by default ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"]
        conditions_stars : list, optional
            Cuts to use to select stars, by default ["extendedness==0"]
        binned_quantity : str, optional
            Quantities to bin the catalog, by default "mag_i_cModel"
        bins : list, optional
            Bins of 'quantities', by default [0, 17, 18, 20, 22, 24]
        """
        self.config_file = config_file
        with open(self.config_file, "r") as f:
            output = yaml.safe_load(f)
        self.name = output.get("name")
        self.openDC2 = call_dc2.OpenDC2(name=self.name)
        self.theta_bins = np.array(output.get("theta_bins"))
        # self.large_theta_bins = np.array(output.get("large_theta_bins"))
        tract_list = output.get("tract_list")
        self.quantities = output.get("quantities")
        self.conditions = output.get("conditions")
        self.conditions_galaxies = output.get("conditions_galaxies")
        self.conditions_stars = output.get("conditions_stars")
        self.binned_quantity = output.get("binned_quantity")
        self.bins = output.get("bins")
        self.conditions_unique_stars = output.get("conditions_unique_stars")
        self.outpath = output.get("outpath")
        if tract_list is None:
            self.tract_list = [
                5074,
                5073,
                5072,
                5071,
                5070,
                5069,
                5068,
                5067,
                5066,
                5065,
                4860,
                4859,
                4858,
                4857,
                4856,
                4855,
                4854,
                4853,
                4852,
                4851,
                4850,
                4648,
                4647,
                4646,
                4645,
                4644,
                4643,
                4642,
                4641,
                4640,
                4639,
                4638,
                4637,
                4636,
                4441,
                4440,
                4439,
                4438,
                4437,
                4436,
                4435,
                4434,
                4433,
                4432,
                4431,
                4430,
                4429,
                4236,
                4235,
                4234,
                4233,
                4232,
                4231,
                4230,
                4229,
                4228,
                4227,
                4226,
                4225,
                4224,
                4035,
                4034,
                4033,
                4032,
                4031,
                4030,
                4029,
                4028,
                4027,
                4026,
                4025,
                4024,
                4023,
                3837,
                3836,
                3835,
                3834,
                3833,
                3832,
                3831,
                3830,
                3829,
                3828,
                3827,
                3826,
                3825,
                3643,
                3642,
                3641,
                3640,
                3639,
                3638,
                3637,
                3636,
                3635,
                3634,
                3633,
                3632,
                3631,
                3453,
                3452,
                3451,
                3450,
                3449,
                3448,
                3447,
                3446,
                3445,
                3444,
                3443,
                3442,
                3441,
                3268,
                3267,
                3266,
                3265,
                3264,
                3263,
                3262,
                3261,
                3260,
                3259,
                3258,
                3257,
                3256,
                3086,
                3085,
                3084,
                3083,
                3082,
                3081,
                3080,
                3079,
                3078,
                3077,
                3076,
                3075,
                3074,
                2908,
                2907,
                2906,
                2905,
                2904,
                2903,
                2902,
                2901,
                2900,
                2899,
                2898,
                2897,
                2896,
            ]
        else:
            self.tract_list = tract_list

    def _open_catalogs(
        self,
        tract=None,
        conditions=None,
        conditions_galaxies=None,
        conditions_stars=None,
        compute_galaxies=True,
        compute_stars=True,
        bins=None,
    ):
        if tract is None:
            tract = self.tract_list
        if bins is None:
            bins = self.bins
        if conditions is None:
            conditions = self.conditions
        if conditions_galaxies is None:
            conditions_galaxies = self.conditions_galaxies
        if conditions_stars is None:
            conditions_stars = self.conditions_stars
        cdt_nan, cdt2 = False, []
        for cdt in conditions_stars:
            if "nan" in cdt:
                cdt_nan = True
                var_nan = cdt.split(".")[-1]
            else:
                cdt2.append(cdt)
        conditions_stars = cdt2
        catalogs = []
        print(
            f"\nOpening catalogs in Critical_radius class. Selected config is bins={bins} over {self.binned_quantity}, conditions={conditions},conditions_galaxies={conditions_galaxies}, conditions_stars={conditions_stars}, compute_galaxies={compute_galaxies},compute_stars={compute_stars}.\n"
        )
        if compute_galaxies:
            galaxies = self.openDC2.galaxies(
                quantities=self.quantities,
                conditions=conditions,
                conditions1=conditions_galaxies,
                tract_list=tract,
            )
            galaxies_neighbour = self.openDC2.galaxies_with_neighbours_tracts(
                quantities=self.quantities,
                conditions=conditions,
                conditions1=conditions_galaxies,
                tract_list=tract,
            )
            print(
                f"\nGalaxies :\n{galaxies[:5]},\n lenght={len(galaxies)};\nGalaxies_neighbour :\n{galaxies_neighbour[:5]},\n lenght={len(galaxies_neighbour)}.\n"
            )
            catalogs.append(galaxies)
            catalogs.append(galaxies_neighbour)
        if compute_stars:
            stars = self.openDC2.stars(
                quantities=self.quantities,
                conditions=conditions,
                conditions1=conditions_stars,
                tract_list=tract,
            )
            if cdt_nan:
                stars = stars[np.isnan(stars[f"{var_nan}"])]
            if len(bins) > 1:
                print(f"There are {len(bins)-1} bins for stars.")
                stars = self.openDC2.bin_cat(
                    stars, quantities=self.binned_quantity, bins=bins
                )
            print(f"\nStars :\n {stars},\n lenght={len(stars)}.\n")
            catalogs.append(stars)
        return catalogs

    def _compute_density_ratio(self, galaxies, galaxies_neighbour, stars):
        density_tracts, number_galaxies, tract_density = (
            np.zeros((len(self.bins) - 1, len(self.theta_bins))),
            0,
            0,
        )
        for i in range(len(stars)):
            compute_density_ = Calculate_density(galaxies, galaxies_neighbour, stars[i])
            density_tracts[i] += compute_density_.get_density_around_stars(
                theta_bins=self.theta_bins
            )
        number_galaxies += len(galaxies["ra"])
        tract_density += compute_density_._tract_density()
        density = density_tracts
        density_dc2 = tract_density
        density_ratio = density / density_dc2
        return density_ratio

    def get_tract_density_ratio(self, tract):
        print(
            f"Now opening catalogs by tracts to compute density_ratio by tract. tract = {tract}."
        )
        galaxies, galaxies_neighbour, stars = self._open_catalogs(tract=tract)
        density_ratio = self._compute_density_ratio(galaxies, galaxies_neighbour, stars)
        return density_ratio

    def get_density_ratio(self):
        if len(self.tract_list) > 1:
            print(
                f"Lenght of tract_list = {len(self.tract_list)} > 1. Performing multiprocessing : see log files for density ratio by tract details."
            )
            mpc = Multiprocessing(
                config_file=self.config_file, tract_list=self.tract_list
            )
            density_ratios = mpc.slurm_submit(outpath=self.outpath)
            density_ratio = np.nanmean(density_ratios, axis=0)
        else:
            print(f"Calculating density_ratio of tract = {self.tract_list}.")
            density_ratio = self.get_tract_density_ratio(self.tract_list)
        return density_ratio

    def get_critical_radius(
        self,
        density_ratio=None,
        critical_density=0.9,
        theta_bins=None,
        binned=None,
        unique=False,
    ):
        """get_critical_radius gets the radius for which to cut for to get a density ratio > critical_density around each bright star

        Parameters
        ----------
        density_ratio : array, optional
            Density_ratio(r) for each bin of stars. if None : computes it using previous function, by default None
        critical_density : float, optional
            Critical density to get around each star, by default 0.9

        Returns
        -------
        array
            radius to cut for each bin of star
        """
        if theta_bins is None:
            theta_bins = self.theta_bins
        if not unique:
            print(f"Calculating critical_radius with theta_bins={theta_bins}.")
        if density_ratio is None:
            density_ratio = self.get_density_ratio()
        if binned is None:
            critical_radius_value = []
            for i in range(len(self.bins) - 1):
                if len(np.where(density_ratio[i] >= critical_density)[0]) > 0:
                    critical_radius_value.append(
                        round(
                            theta_bins[
                                np.where(density_ratio[i] >= critical_density)[0][0]
                            ],
                            2,
                        )
                    )
                else:
                    critical_radius_value.append(np.nan)
        else:
            if unique == False:
                critical_radius_value = round(
                    theta_bins[np.where(density_ratio >= critical_density)[0][0]], 2
                )
            else:
                if len(theta_bins[np.where(density_ratio >= critical_density)[0]]) > 1:
                    critical_radius_value = theta_bins[
                        np.where(density_ratio >= critical_density)[0][0]
                    ]
                else:
                    print(
                        f"Issue detected, setting critical radius value to 0. Some density_ratios might be < {critical_density} for all thetas."
                    )
                    critical_radius_value = 0

        return critical_radius_value

    def get_unique_tract_density_ratio(self, tract):
        print(f"\nConditions in unique stars : see called method bellow :\n")
        galaxies, galaxies_neighbour, stars = self._open_catalogs(
            tract=tract,
            conditions=None,
            conditions_galaxies=None,
            conditions_stars=None,
            bins=None,
        )
        density_brighter = []
        ra, dec = [], []
        if len(self.bins) > 1:
            for i in range(len(self.bins) - 1):
                stars_bin = stars[i]
                density_brighter_bin = np.zeros((len(stars_bin), len(self.theta_bins)))
                compute_density_ = Calculate_density(
                    galaxies, galaxies_neighbour, stars_bin
                )
                density_brighter_bin += compute_density_.get_density_around_unique_star(
                    theta_bins=self.theta_bins
                )
                tract_density = compute_density_._tract_density()
                density = density_brighter_bin
                density_dc2 = tract_density
                density_ratio = density / density_dc2
                density_brighter.append(density_ratio)
                ra.append(stars_bin["ra"])
                dec.append(stars_bin["dec"])
        else:
            stars_bin = stars
            density_brighter_bin = np.zeros((len(stars_bin), len(self.theta_bins)))
            compute_density_ = Calculate_density(
                galaxies, galaxies_neighbour, stars_bin
            )
            density_brighter_bin += compute_density_.get_density_around_unique_star(
                theta_bins=self.theta_bins
            )
            tract_density = compute_density_._tract_density()
            density = density_brighter_bin
            density_dc2 = tract_density
            density_ratio = density / density_dc2
            density_brighter.append(density_ratio)
            ra.append(stars_bin["ra"])
            dec.append(stars_bin["dec"])
        return ra, dec, density_brighter

    def get_unique_density_ratio(self):
        if len(self.tract_list) > 1:
            print(
                f"Lenght of tract_list ={len(self.tract_list)} > 1. Performing multiprocessing : see log files for density ratio by tract details."
            )
            mpc = Multiprocessing(tract_list=self.tract_list, unique=True)
            density_ratios_table = mpc.slurm_submit(outpath=self.outpath)
            ra, dec, density_ratio = (
                density_ratios_table["ra"],
                density_ratios_table["dec"],
                density_ratios_table["density_ratio"],
            )
        else:
            print(f"Calculating density_ratio of tract={self.tract_list}.")
            ra, dec, density_ratio = self.get_unique_tract_density_ratio(
                self.tract_list
            )
        return ra, dec, density_ratio

    def get_unique_critical_radius(
        self,
        ra=None,
        dec=None,
        density_ratio=None,
        critical_density=0.9,
        theta_bins=None,
    ):
        if theta_bins is None:
            theta_bins = self.theta_bins
        if density_ratio is None:
            ra, dec, density_ratio = self.get_unique_density_ratio()
        critical_radius_value = np.zeros(len(density_ratio))
        for i in range(len(critical_radius_value)):
            unique_critical_radius = self.get_critical_radius(
                density_ratio=density_ratio[i],
                theta_bins=theta_bins,
                critical_density=critical_density,
                unique=True,
                binned=True,
            )
            critical_radius_value[i] = unique_critical_radius
        return ra, dec, critical_radius_value
