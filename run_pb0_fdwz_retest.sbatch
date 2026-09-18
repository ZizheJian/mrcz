#!/bin/bash
#SBATCH --job-name=fdwz-retest
#SBATCH --account=cis240669
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=00:05:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
srun --ntasks=1 --cpus-per-task=64 --cpu-bind=cores python3 -u py_funcs/benchmark_fdwz_selector.py --dataset pb0 --fields P --ebs 1e-4,1e-5 --max-cpus 64 --repeats 3 --sample-interval 128 --strategies linear-proxy,SZ3-4D --output results/RD/pb0_fdwz_retest_clean_20260905.csv
