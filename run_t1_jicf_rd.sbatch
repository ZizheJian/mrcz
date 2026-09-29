#!/bin/bash
#SBATCH --job-name=mrcz-jicf-rd
#SBATCH --account=cis240669
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=235G
#SBATCH --time=01:10:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export MRCZ_T1_DATASET=jicf0
export MRCZ_T1_MAX_CPUS=128
export MRCZ_T1_BUILD_BEFORE_TEST=false
export MRCZ_T1_METHODS=FDWZ
export MRCZ_T1_FIELDS=UX,UY,UZ,P,S00
export MRCZ_T1_EBS=1e-4,1e-5,1e-6,1e-7,1e-8
export MRCZ_T1_FDWZ_SAMPLE_INTERVAL=128
TIMEFORMAT=$'\nTotal wall time: %R seconds\nTotal CPU time: %U user + %S system seconds'
time srun --ntasks=1 --cpus-per-task=128 --cpu-bind=cores python3 -u t1_start_c_code.py
python3 py_funcs/merge_method_csv.py /anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/viral-flamespeed results/RD/jicf0 --method FDWZ --pattern 'jicf0.f00000_*.csv'
python3 py_funcs/plot_large_rd.py jicf0
python3 py_funcs/plot_large_throughput.py jicf0 --x-axis eb
python3 py_funcs/plot_large_throughput.py jicf0
python3 py_funcs/summarize_rd.py
