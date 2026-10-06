#!/bin/bash
#SBATCH --job-name=mrcz-vis-tune
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --time=00:15:00
#SBATCH --output=results/visualization/slurm-%x-%j.out
#SBATCH --error=results/visualization/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
mkdir -p results/visualization
if [[ ! -e results/visualization/pb0.f00001 && ! -L results/visualization/pb0.f00001 ]];then
    ln -s /anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/hyuan/pb0.f00001 results/visualization/pb0.f00001
fi
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
config_name="${CONFIG_NAME:-visualization_cr_tuning.json}"
srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores python3 -u t1_start_c_code.py "${config_name}"
