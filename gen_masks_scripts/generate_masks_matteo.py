import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import hpgeom as hpg
import healsparse as hsp
import healpy as hp
from astropy.io import fits
import yaml

config_file = "/pbs/home/n/namourou/workspace/side_codes/clusters/dc2/desc_stellar_masks_package/config/config.yaml"
with open(config_file, "r") as f:
    output = yaml.safe_load(f)
outpath = output.get("outpath")

print("Now convert to healpixels")
nside_hp = int(4096 * 2)
m = hsp.HealSparseMap.read(outpath + "bo_masks.hs")
# methods = ['mean', 'median', 'std', 'max', 'min', 'sum', 'prod']
method = "sum"
# for method in methods:
hp_map = m.generate_healpix_map(nside=nside_hp, reduction="sum")
print("Healpix map generated")
ra_hp, dec_hp = hp.pix2ang(nside_hp, np.arange(len(hp_map)), lonlat=True, nest=True)
print("Positions found")
ra_max, ra_min, dec_max, dec_min = 72, 52, -26, -46
ra_hp_s = ra_hp[
    (ra_hp < ra_max) & (ra_hp > ra_min) & (dec_hp < dec_max) & (dec_hp > dec_min)
]
dec_hp_s = dec_hp[
    (ra_hp < ra_max) & (ra_hp > ra_min) & (dec_hp < dec_max) & (dec_hp > dec_min)
]
hp_map_s = hp_map[
    (ra_hp < ra_max) & (ra_hp > ra_min) & (dec_hp < dec_max) & (dec_hp > dec_min)
]
ra_hp, dec_hp, hp_map = ra_hp_s, dec_hp_s, hp_map_s

# hp_map = Table({"mask" : hp_map})
# hp_map["mask_bool"] = 1
# hp_map["mask_bool"][hp_map["mask"]/256 > 0.13] = 0
print("Convert to fits standards...")
nside_ratio = hp.nside2npix(131072) / hp.nside2npix(nside_hp)
hp_map_t = Table({"WEIGHT": hp_map[hp_map / nside_ratio > 0.13] / nside_ratio})
# hp_map_t["BOOL"] = 1
hp_map_t["PIXEL"] = hp.ang2pix(
    nside_hp,
    ra_hp[hp_map / nside_ratio > 0.13],
    dec_hp[hp_map / nside_ratio > 0.13],
    lonlat=True,
    nest=True,
)
# hp_map_t["RA"] = ra_hp[hp_map/nside_ratio>0.13]
# hp_map_t["DEC"] = dec_hp[hp_map/nside_ratio>0.13]
header = fits.Header()
# header['RANGE_RA'] = f'{[min(ra_hp), max(ra_hp)]}'
# header['RANGE_DEC'] = f'{[min(dec_hp), max(dec_hp)]}'
header["NSIDE"] = f"{nside_hp}"
header["ORDERING"] = True
header["DOWNGRADE_METHOD"] = "sum for healsparse + > 13% of pixel masked"
print("Now saving")
hdu = fits.BinTableHDU(hp_map_t, header=header)

# Save the FITS file
# hdu.writeto(outpath + f"first_full_mask/bo_masks_healpixels.fits", overwrite=True)
hdu.writeto(outpath + f"bo_masks_matteo.fits", overwrite=True)

# with open(outpath + f"first_full_mask/bo_masks_healpixels.fits", 'wb') as file:
#    pickle.dump(healpix_data, file)
"""
print("Plotting...")
plt.scatter(ra_hp_s, dec_hp_s, c = hp_map_s,s=1)
plt.xlabel("ra", fontsize = 13)
plt.ylabel("dec", fontsize = 13)
plt.savefig(f"/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/plots/study/cat_gen_hp_{tract}_{method}.png")"""
"""fig, ax = plt.subplots(figsize=(8, 6))
sp = skyproj.McBrydeSkyproj(ax=ax)
sp.draw_hspmap(m, valid_mask=True, lon_range=[59.5,59.75], lat_range=[-37.2,-37])
sp.draw_inset_colorbar()
sp.draw_hpxmap(healpix_data["mask"]/256, nest=True)
sp.draw_inset_colorbar()
plt.savefig("/sps/lsst/users/namourou/web/clusters/DC2/bright_objects_masks/plots/study/cat_gen_hp_hsp_{tract}_{method}.png")"""
