import sys

sys.path.append("../bright_objects_masks")
import generate_masks
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
mask_gen = generate_masks.Masks(
    config_file="/pbs/home/n/namourou/workspace/side_codes/clusters/dc2/desc_stellar_masks_package/config/config.yaml"
)
mask = mask_gen.create_healsparse_masks()
mask_gen.write_heaslparse_mask(mask)
