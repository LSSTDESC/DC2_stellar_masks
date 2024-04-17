#!/bin/sh
#sbatch -t  1:00:00 -n 2 --mem 100G
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/namourou/conda_envs/conda_clone_021023/desc_v0/
#python produce_masks.py
#python generate_masks_matteo.py
#echo "env loaded now producing masks."
#python stars_masks.py
#python vizualisation.py
python radius_study_brightest.py