#!/bin/bash
#SBATCH --job-name=mrcz-pb0-rd
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=06:00:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
unset MRCZ_PROFILE_HUFFMAN MRCZ_PROFILE_PHASES
export MRCZ_T1_DATASET=pb0
export MRCZ_T1_BUILD_BEFORE_TEST=false
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time python3 -u t1_start_c_code.py
