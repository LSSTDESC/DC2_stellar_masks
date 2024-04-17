import sys

sys.path.append("../bright_objects_masks")
import radius_study
import generate_masks
import warnings
import healsparse as hsp
from astropy.table import Table
import matplotlib.pyplot as plt
import skyproj
import yaml

warnings.filterwarnings("ignore")

config_file = "/pbs/home/n/namourou/workspace/side_codes/clusters/dc2/desc_stellar_masks_package/config/config.yaml"
with open(config_file, "r") as f:
    output = yaml.safe_load(f)
outpath = output.get("outpath")

mask = hsp.HealSparseMap.read(outpath + "bo_masks.hs")

fig, ax = plt.subplots(figsize=(8, 5))
sp = skyproj.McBrydeSkyproj(ax=ax)
im, lon_raster, lat_raster, values_raster = sp.draw_hspmap(
    mask, vmin=0, vmax=1, lon_range=[52, 72], lat_range=[-45, -25], cmap="Dark2"
)
# sp.plot(ra, dec, "b.", markersize =1) #If you want to vizualize some footprint over masks (galaxies, stars, ...)
plt.annotate(
    "Masks vs stars",
    xy=(0.5, 1.03),
    xycoords="axes fraction",
    xytext=(0, 10),
    textcoords="offset points",
    ha="center",
    va="bottom",
    fontsize=15,
)
plt.colorbar()
plt.savefig(outpath + "masks_map.png")
