#!/bin/sh
sbatch -t  2:00:00 --cpus-per-task=1 -n 2 --mem 200G --partition lsst,htc produce_slurm.sh
