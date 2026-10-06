#!/bin/bash
#SBATCH --job-name=mrcz-weak-jicf-16cpu
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --hint=nomultithread
#SBATCH --time=00:40:00
#SBATCH --output=results/weak_scaling/slurm-%x-%j.out
#SBATCH --error=results/weak_scaling/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
srun --ntasks=1 --cpus-per-task=16 --cpu-bind=cores python3 -u t5_start_weak_scaling.py jicf "$@"
