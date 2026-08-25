#!/bin/bash
#SBATCH --job-name=mrcz-e11
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --output=experiments/sz3_3d_gll/slurm-%x-%j.out
#SBATCH --error=experiments/sz3_3d_gll/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time python3 -u experiments/sz3_3d_gll/e11_compare_sz3_prediction_order.py
