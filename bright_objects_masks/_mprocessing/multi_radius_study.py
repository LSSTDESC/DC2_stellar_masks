import sys

sys.path.append("../")
import radius_study
import numpy as np
import pickle

tract, config_file, outpath, unique = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
print(
    f"\nProcessing tract = {tract},\nWith config file at {config_file},\nSaving at {outpath}\nProcessing by stars : {unique}.\n"
)

critical_radius = radius_study.Critical_radius(config_file=config_file)
if eval(unique):
    print("\nInitialization of processing by star.\n")
    ra, dec, density_ratio = critical_radius.get_unique_tract_density_ratio(tract=tract)
    results = {"ra": ra, "dec": dec, "density_ratio": density_ratio}
    with open(outpath + f"{tract}_density_ratio_unique.pickle", "wb") as handle:
        pickle.dump(results, handle)
else:
    print("\nInitialization of processing by bins.\n")
    density_ratio = critical_radius.get_tract_density_ratio(tract=tract)
    print(f"\nResult of processing for tract {tract} : \n{density_ratio}.\n")
    np.savetxt(outpath + f"{tract}_density_ratio.txt", density_ratio)
