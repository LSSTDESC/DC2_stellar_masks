import sys

sys.path.append("../../bright_objects_masks")
import radius_study
import warnings

warnings.filterwarnings("ignore")

print("\nFirst basic instantiation.\n")
critical_radius = radius_study.Critical_radius(
    config_file="/pbs/home/n/namourou/workspace/side_codes/clusters/dc2/desc_stellar_masks_package/examples/test_package/config_unique_stars_multiple_tracts.yaml"
)
print("\nGet density ratios.\n")
ra, dec, density_ratio = critical_radius.get_unique_density_ratio()
print(ra, dec, density_ratio)
print("\nGet critical radiuses.\n")
ra2, dec2, crit_radius = critical_radius.get_unique_critical_radius(
    ra=ra, dec=dec, density_ratio=density_ratio
)
print(crit_radius)
print("\nOk for nominal processing.")
