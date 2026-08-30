#!/bin/bash
#SBATCH --job-name=mrcz-weak-jicf
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread
#SBATCH --time=06:00:00
#SBATCH --output=experiments/mrcz_scaling/slurm-%x-%j.out
#SBATCH --error=experiments/mrcz_scaling/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export MRCZ_DISABLE_CONSTANT_ENCODING=1
unset OMP_PLACES OMP_PROC_BIND MRCZ_PROFILE_MEMORY MRCZ_PROFILE_HUFFMAN MRCZ_PROFILE_PHASES || true
srun --ntasks=1 --cpus-per-task=128 --cpu-bind=cores python3 -u experiments/mrcz_scaling/run_weak_scaling.py jicf --append-existing --repair-log experiments/mrcz_scaling/slurm-mrcz-weak-jicf-20205857.out --methods SZ3-4D-uniform,SZ2-3D-uniform-fused,QoZ-4D-uniform,ZFP-3D-uniform-fused,SPERR-2D-uniform
