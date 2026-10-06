#!/bin/bash
#SBATCH --job-name=mrcz-small-jicf-split
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --hint=nomultithread
#SBATCH --time=00:45:00
#SBATCH --output=results/small_core_throughput/slurm-%x-%j.out
#SBATCH --error=results/small_core_throughput/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
mkdir -p results/small_core_throughput
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
config_name="${CONFIG_NAME:?CONFIG_NAME must name one method-specific JICF configuration}"
srun --ntasks=1 --cpus-per-task=48 --cpu-bind=cores python3 -u t1_start_c_code.py "${config_name}" --skip-build
