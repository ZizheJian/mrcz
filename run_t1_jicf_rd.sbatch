#!/bin/bash
#SBATCH --job-name=mrcz-jicf-rd
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=10:00:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
unset MRCZ_PROFILE_HUFFMAN
export MRCZ_T1_DATASET=jicf0
export MRCZ_T1_MAX_CPUS=128
export MRCZ_T1_BUILD_BEFORE_TEST=false
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time python3 -u t1_start_c_code.py
cp /anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000_*.csv results/RD/jicf0/
