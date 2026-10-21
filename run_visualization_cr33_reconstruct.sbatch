#!/bin/bash
#SBATCH --job-name=mrcz-vis-T-cr33-recon
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --hint=nomultithread
#SBATCH --time=00:10:00
#SBATCH --output=results/visualization/slurm-%x-%j.out
#SBATCH --error=results/visualization/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
mkdir -p results/visualization/cr33_reconstruction
mkdir -p /anvil/projects/x-cis240192/x-zjian1-in-cis240192/mrcz-visualization/cr33
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores python3 -u t1_start_c_code.py visualization_cr33_reconstruct_T.json --skip-build
