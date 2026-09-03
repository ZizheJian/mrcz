#!/bin/bash
#SBATCH --job-name=mrcz-pb0-rd
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
unset MRCZ_PROFILE_HUFFMAN
export MRCZ_T1_DATASET=pb0
export MRCZ_T1_MAX_CPUS=64
export MRCZ_T1_BUILD_BEFORE_TEST=false
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time srun --ntasks=1 --cpus-per-task=64 --cpu-bind=cores python3 -u t1_start_c_code.py
cp /anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001_*.csv results/RD/pb0/
