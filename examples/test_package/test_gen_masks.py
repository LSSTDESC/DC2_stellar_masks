import sys

sys.path.append("../../bright_objects_masks")
import generate_masks
import warnings

warnings.filterwarnings("ignore")

print("\nFirst basic instantiation.\n")
mask_gen = generate_masks.Masks(
    config_file="/pbs/home/n/namourou/workspace/side_codes/clusters/dc2/desc_stellar_masks_package/examples/test_package/config_multiple_tracts.yaml"
)  # , density_ratio = density_glob)
mask = mask_gen.create_healsparse_masks()
mask_gen.write_heaslparse_mask(mask)
print("\nOk for nominal processing.")
