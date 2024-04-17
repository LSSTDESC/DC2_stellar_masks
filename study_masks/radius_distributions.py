import sys

sys.path.append("../bright_objects_masks")
import radius_study
import call_dc2
import warnings
import numpy as np
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

critical_radius = radius_study.Critical_radius(
    name="dc2_object_run2.2i_dr6_v2_with_addons_v2",
    theta_bins=np.logspace(np.log10(0.5), np.log10(500), 200),
    tract_list=None,
    quantities=["ra", "dec", "mag_i_cModel", "mag_i_truth"],
    conditions=["mag_i_truth<27"],
    conditions_galaxies=[
        "extendedness==1",
        "mag_i_cModel>17",
        "mag_i_cModel<25.3",
        "clean==True",
    ],
    conditions_stars=["extendedness!=1.0"],
    binned_quantity="mag_i_cModel",
    bins=[15, 16, 17],
)
mag_bins = [15, 16, 17]
density_glob = critical_radius.get_density_ratio()
print(density_glob)
np.savetxt(
    "/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/study_masks/density_ratio_truemag.txt",
    density_glob,
)

for i in range(len(density_glob)):
    plt.figure()
    plt.plot(
        np.logspace(np.log10(0.5), np.log10(500), 200),
        density_glob[i],
        linestyle="",
        marker="+",
    )
    plt.xscale("log")
    plt.savefig(
        f"/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/study_masks/density_vs_theta{i}.png"
    )
    plt.close

print("now processing unique bright stars")
for i in range(len(mag_bins) - 1):
    ra, dec, density_unique_bin = critical_radius.get_unique_tract_density_ratio(
        "3830",
        conditions=[""],
        conditions_stars=[
            "extendedness!=1",
            "clean==False",
            "nan.mag_i_cModel",
            f"mag_i_truth<{mag_bins[i+1]}",
            f"mag_i_truth>{mag_bins[i]}",
        ],
        larger_theta_bins=np.logspace(np.log10(0.5), np.log10(500), 200),
    )
    print(density_unique_bin)
    np.savetxt(
        f"/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/study_masks/density_ratio_truemag_bin{i}.txt"
    )
    radius = critical_radius.get_unique_critical_radius(
        ra,
        dec,
        density_unique_bin,
        larger_theta_bins=np.logspace(np.log10(0.5), np.log10(500), 200),
    )
    plt.figure()
    plt.hist(radius, bins=25)
    plt.savefig(
        f"/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/study_masks/radius_distrib_truemag{i}.png"
    )
    plt.close()
