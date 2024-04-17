import sys

sys.path.append("bright_objects_masks")
import radius_study
import call_dc2
import warnings
import numpy as np
from astropy.table import Table

warnings.filterwarnings("ignore")


critical_radius = radius_study.Critical_radius(tract_list=None)  # ,4438])

ra, dec, brighter_critical_density_value = critical_radius.get_unique_density_ratio()


density_ratio = Table(
    {"ra": ra, "dec": dec, "density_ratio": brighter_critical_density_value}
)
density_ratio.write(
    "/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/bo_custom_ratio_multiple.fits",
    overwrite=True,
)


ra, dec, brighter_critical_radius_value = critical_radius.get_unique_critical_radius(
    ra=ra, dec=dec, density_ratio=brighter_critical_density_value
)
print(brighter_critical_radius_value)
density_ratio["radius"] = brighter_critical_radius_value
density_ratio.write(
    "/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/bo_custom_ratio_radius_multiple.fits",
    overwrite=True,
)
