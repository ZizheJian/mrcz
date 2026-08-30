#!/bin/bash
#SBATCH --job-name=mrcz-scaling
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread
#SBATCH --time=00:30:00
#SBATCH --output=experiments/mrcz_scaling/slurm-%x-%j.out
#SBATCH --error=experiments/mrcz_scaling/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
unset OMP_PLACES || true
unset OMP_PROC_BIND || true
unset MRCZ_PROFILE_PHASES || true

cmake -S c_codes -B build -DCMAKE_INSTALL_PREFIX="${SLURM_SUBMIT_DIR}/install"
cmake --build build --parallel 8
cmake --install build
hostname
lscpu
sha256sum install/bin/mrcz_compress install/bin/mrcz_decompress
srun --ntasks=1 --cpus-per-task=128 --cpu-bind=cores python3 -c "import os;print('Slurm step CPU affinity=',len(os.sched_getaffinity(0)),flush=True)"
srun --ntasks=1 --cpus-per-task=128 --cpu-bind=cores python3 -u experiments/mrcz_scaling/run_scaling.py
