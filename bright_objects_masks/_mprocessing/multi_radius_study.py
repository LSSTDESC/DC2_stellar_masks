import sys
sys.path.append('../')
import radius_study
import yaml
import numpy as np
import os

outpath = os.path.dirname(__file__) + "/logs/"

tract, config_file = sys.argv[1], sys.argv[2]

with open(config_file,'r') as f:
    output = yaml.safe_load(f)
name = output.get('name')
theta_bins = np.array(output.get('theta_bins'))
quantities = output.get('quantities')
conditions = output.get('conditions')
conditions_galaxies = output.get('conditions_galaxies')
conditions_stars = output.get('conditions_stars')
binned_quantity = output.get('binned_quantity')
bins = output.get('bins')

critical_radius = radius_study.Critical_radius(name=name, theta_bins=theta_bins,tract_list=tract, quantities=quantities, conditions=conditions,
                                                conditions_galaxies=conditions_galaxies, conditions_stars=conditions_stars, binned_quantity=binned_quantity, bins=bins)
density_ratio = critical_radius.get_tract_density_ratio(tract=tract)
np.savetxt(outpath + f"{tract}_density_ratio.txt", density_ratio)
