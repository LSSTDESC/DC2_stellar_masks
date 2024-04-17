#!/bin/sh
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/namourou/conda_envs/conda_clone_021023/desc_v0/
python produce_masks.py
python generate_masks_matteo.py
python vizualisation.py