import yaml
import numpy as np

parameters = {
    "name": "dc2_object_run2.2i_dr6_v2_with_addons_v2",
    "theta_bins": np.logspace(np.log10(0.5), np.log10(50), 50).tolist(),
    "tract_list": [3830, 4433],  # Full DC2 catalog
    "quantities": ["ra", "dec", "mag_i_cModel"],
    "conditions": None,  #'Basic quality cuts'
    "conditions_galaxies": ["extendedness==1", "mag_i_cModel>17", "mag_i_cModel<25.3"],
    "conditions_stars": ["extendedness==0"],
    "binned_quantity": "mag_i_cModel",
    "bins": [0, 17, 18, 20, 22, 24],
    "density_ratio": None,  # The function will auto compoute it
    "critical_density": 0.9,
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
