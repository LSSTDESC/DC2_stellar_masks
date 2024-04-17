import time
import subprocess
import os
import sys

sys.path.append("../")
import numpy as np
import yaml
import pickle

current_dir = os.path.dirname(__file__)


class Multiprocessing:
    def __init__(
        self,
        config_file=f"{current_dir}/../../config/config.yaml",
        tract_list=None,
        unique=False,
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
        self.theta_bins = np.array(output.get("theta_bins"))
        self.bins = output.get("bins")
        self.tract_list = tract_list
        self.unique = unique

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
        slurm_mem = 149
        res, job_list, slurm_output_path_list = [], [], []
        if outpath is None:
            outpath = current_dir + "/logs/"
        self._mkdir(outpath)
        print("\nNow submitting jobs for multiprocessing.\n")
        for tract in self.tract_list:
            slurm_output_path = f"{outpath}{func_name}_{tract}.out"
            slurm_output_path_list.append(slurm_output_path)
            cmd = f"sbatch --job-name=submit_func_{func_name}_{tract} -t  0-15:00 -n 2 --mem {slurm_mem}G --partition lsst,htc -D {current_dir} -L sps -o {slurm_output_path} <<?\n"
            cmd += "#!/usr/bin/bash\n"
            cmd += f"source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh\n"
            cmd += f"conda activate /sps/lsst/users/namourou/conda_envs/conda_clone_021023/desc_v0/\n"
            cmd += f"python {current_dir}/{func_name}.py {tract} {self.config_file} {outpath} {self.unique}"
            res = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, check=True
            )
            job_id = str(res.stdout).split("batch job")[1].split("\\")[0]
            job_list.append(job_id.split("\n")[0])
        if not self.unique:
            print("\nUsing nominal processing.\n")
            result = np.zeros(
                (len(self.tract_list), len(self.bins) - 1, len(self.theta_bins))
            )  # create empty array to store results
        for i, tract in enumerate(self.tract_list):
            self._check_slurm_job(job_list[i])
            if self.unique:
                print("\nProcessing by star.\n")
                if i == 0:
                    with open(
                        outpath + f"{tract}_density_ratio_unique.pickle", "rb"
                    ) as f:
                        result = pickle.load(f)
                else:
                    with open(
                        outpath + f"{tract}_density_ratio_unique.pickle", "rb"
                    ) as f:
                        curr_result = pickle.load(f)
                    result["ra"] = result["ra"] + curr_result["ra"]
                    result["dec"] = result["dec"] + curr_result["dec"]
                    result["density_ratio"] = (
                        result["density_ratio"] + curr_result["density_ratio"]
                    )
            else:
                result[i] = np.loadtxt(outpath + f"{tract}_density_ratio.txt")
        cmd = f"cp {self.config_file} {outpath}\n"
        res = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, check=True
        )
        return result
