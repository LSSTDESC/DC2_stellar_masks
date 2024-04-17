import yaml
import numpy as np

parameters = {
    "name": "dc2_object_run2.2i_dr6_v2_with_addons_v2",
    "theta_bins": np.logspace(np.log10(0.5), np.log10(500), 200).tolist(),
    "tract_list": None,  # Full DC2 catalog set to None or use list even if one tract ex: ["3830"]
    "quantities": ["ra", "dec", "mag_i_cModel", "mag_i_truth"],
    "conditions": ["mag_i_truth<27"],  #'Basic quality cuts'
    "conditions_galaxies": ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"],
    "conditions_stars": ["extendedness!=1", "clean==False", "nan.mag_i_cModel"],
    "binned_quantity": "mag_i_truth",
    "bins": [0, 10, 11, 12, 13, 14, 15, 16],
    "density_ratio": None,  # The function will compute it for us
    "critical_density": 0.9,
    "outpath": "/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/catalogs/full_dc2_masks_mag_i_truth/",
    "unique_stars": False,
    "nside_coverage": 32,
    "nside_sparse": 131072,
}
yaml_output = yaml.dump(parameters, sort_keys=False)
print(yaml_output)
with open(
    "config.yaml",
    "w",
) as f:
    yaml.dump(parameters, f, sort_keys=False)
    print("Written to file successfully")
