#!/bin/bash
#SBATCH --job-name=mrcz-t1
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time srun --ntasks=1 --cpus-per-task=128 --cpu-bind=cores python3 -u t1_start_c_code.py
