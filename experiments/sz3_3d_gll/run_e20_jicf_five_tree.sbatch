#!/bin/bash
#SBATCH --job-name=mrcz-e20b
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:45:00
#SBATCH --output=experiments/sz3_3d_gll/slurm-%x-%j.out
#SBATCH --error=experiments/sz3_3d_gll/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export MRCZ_T1_BUILD_BEFORE_TEST=false
export MRCZ_T1_DATASET=jicf0
export MRCZ_T1_METHODS=SZ3-3D-GLL-experimental
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time python3 -u t1_start_c_code.py
