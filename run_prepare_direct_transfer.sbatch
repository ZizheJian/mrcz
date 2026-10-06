#!/bin/bash
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --job-name=mrcz-direct-data
#SBATCH --output=results/end_to_end/slurm-%x-%j.out
#SBATCH --error=results/end_to_end/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
"${SLURM_SUBMIT_DIR}/mrcz_env/bin/python" py_funcs/prepare_direct_transfer.py
