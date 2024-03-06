import time
import subprocess
import os
import sys

sys.path.append("../")
import call_dc2
import numpy as np
import yaml

current_dir = os.path.dirname(__file__)


class Multiprocessing:
    def __init__(
        self, config_file=f"{current_dir}/../../config/config.yaml", tract_list=None
    ):
        """__init__ class which submit slurm job to produce masks (in particular for density ratio computation)

        Parameters
        ----------
        config_file : str, optional
            Path of config_file, by default f'{current_dir}/../../config/config.yaml'.
        tract_list : list, optional
            List of tracts to parallelize, by default None = full catalog.
        """
        self.config_file = config_file
        with open(self.config_file, "r") as f:
            output = yaml.safe_load(f)
        self.name = output.get("name")
        self.openDC2 = call_dc2.OpenDC2(name=self.name)
        self.theta_bins = np.array(output.get("theta_bins"))
        self.tract_list = tract_list
        self.quantities = output.get("quantities")
        self.conditions = output.get("conditions")
        self.conditions_galaxies = output.get("conditions_galaxies")
        self.conditions_stars = output.get("conditions_stars")
        self.binned_quantity = output.get("binned_quantity")
        self.bins = output.get("bins")

    def _mkdir(self, dir):
        if not os.path.exists(dir):
            os.mkdir(dir)

    def _check_slurm_job(self, job_id):
        """check_slurm_job Checks if submitted jobs are finished.

        Parameters
        ----------
        job_id : str
            ID of the slurm job.
        """
        while True:
            res = subprocess.run(
                ["squeue", "-j", job_id], capture_output=True, text=True
            )
            output = res.stdout
            if job_id not in output:
                break
            time.sleep(30)

    def slurm_submit(
        self, func_name="multi_radius_study", outpath=current_dir + "/logs/"
    ):
        slurm_mem = 64
        res, job_list, slurm_output_path_list = [], [], []
        self._mkdir(outpath)
        for tract in self.tract_list:
            slurm_output_path = f"{outpath}{func_name}_{tract}.out"
            slurm_output_path_list.append(slurm_output_path)
            cmd = f"sbatch --job-name=submit_func_{func_name}_{tract} -t  0-15:00 -n 2 --mem {slurm_mem}G -D {current_dir} -L sps -o {slurm_output_path} <<?\n"
            cmd += "#!/usr/bin/bash\n"
            cmd += f"source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh\n"
            cmd += f"conda activate /sps/lsst/users/namourou/conda_envs/conda_clone_021023/desc_v0/\n"
            cmd += f"python {current_dir}/{func_name}.py {tract} {self.config_file}"
            res = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, check=True
            )
            job_id = str(res.stdout).split("batch job ")[1].split("\\")[0]
            job_list.append(job_id.split("\n")[0])
        result = np.zeros(
            (len(self.tract_list), len(self.bins) - 1, len(self.theta_bins))
        )  # create empty arry to store results
        for i, tract in enumerate(self.tract_list):
            self._check_slurm_job(job_list[i])
            result[i] = np.loadtxt(outpath + f"{tract}_density_ratio.txt")
        return result
