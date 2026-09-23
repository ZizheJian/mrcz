#!/bin/bash
#SBATCH --job-name=prediction-entropy
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --time=01:00:00
#SBATCH --output=results/prediction_locality/slurm-%x-%j.out
#SBATCH --error=results/prediction_locality/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
srun --ntasks=1 --cpus-per-task=8 --cpu-bind=cores python3 -u t1_start_c_code.py pb0_fdwz_prediction_locality
python3 -u t3_test_cubic_face_prediction.py
