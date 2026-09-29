#!/bin/bash
#SBATCH --job-name=fdwz-tree-qb
#SBATCH --account=cis240669
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=235G
#SBATCH --time=00:05:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
srun --ntasks=1 --cpus-per-task=64 --cpu-bind=cores install/bin/mrcz_compress -format nekrs -i /anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/hyuan/pb0.f00001 -M REL 0.0001 -method FDWZ-ablation -max-cpus 64 -partition-num 64 -fields '[P]'
