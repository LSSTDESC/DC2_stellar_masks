import sys

sys.path.append("../../bright_objects_masks")
import radius_study
import warnings

warnings.filterwarnings("ignore")

print("\nFirst basic instantiation.\n")
critical_radius = radius_study.Critical_radius(config_file="config_one_tract.yaml")
print("\n_open_catalogs method\n")
stars1, galaxies1, galaxies_w_neighbour1 = critical_radius._open_catalogs()
print("\nget_density_ratio method with single tract.\n")
density_ratio = critical_radius.get_density_ratio()
print("\nget critical radius.\n")
crit_radius = critical_radius.get_critical_radius(density_ratio=density_ratio)
print(crit_radius)
print("\nNow trying multiprocessing with several tracts\n")
print("\nSecond instantiation\n")
critical_radius_m = radius_study.Critical_radius(
    config_file="/pbs/home/n/namourou/workspace/side_codes/clusters/dc2/desc_stellar_masks_package/examples/test_package/config_multiple_tracts.yaml"
)
print("\nGet critical radiuses\n")
crit_radius_m = critical_radius_m.get_critical_radius()
print(crit_radius_m)
print("\nOk for nominal processing.")
