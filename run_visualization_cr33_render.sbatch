#!/bin/bash
#SBATCH --job-name=mrcz-vis-T-cr33-render
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
output_dir="results/visualization/cr33_slices"
raw_dir="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/mrcz-visualization/cr33"
mkdir -p "${output_dir}"
module load paraview/5.10.1
export PYTHONPATH="${SLURM_SUBMIT_DIR}/paraview_env${PYTHONPATH:+:${PYTHONPATH}}"
export PB0_FIELD=T
unset PB0_FIELD_VALUES_FILE
export PB0_SELECTED_ELEMENTS_OUTPUT="${output_dir}/pb0_T_x0_selected.npy"
export PB0_X0_SLICE_DATA="${output_dir}/pb0_T_x0_original.npz"
srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores pvbatch py_funcs/render_pb0_x0_subset.py
unset PB0_SELECTED_ELEMENTS_OUTPUT
export PB0_SELECTED_ELEMENTS_INPUT="${output_dir}/pb0_T_x0_selected.npy"
methods=(MRCZ FDWZ SZ3-4D SZ2-3D-fused QoZ-4D ZFP-3D-fused SPERR-2D)
for method in "${methods[@]}";do
    raw_path="${raw_dir}/pb0_T_${method}.raw"
    if [[ ! -f "${raw_path}" ]];then
        echo "Missing reconstructed field: ${raw_path}" >&2
        exit 1
    fi
    export PB0_FIELD_VALUES_FILE="${raw_path}"
    export PB0_X0_SLICE_DATA="${output_dir}/pb0_T_x0_${method}.npz"
    srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores pvbatch py_funcs/render_pb0_x0_subset.py
    rm -f "${raw_path}"
done
module purge
unset PYTHONPATH PB0_FIELD PB0_FIELD_VALUES_FILE PB0_SELECTED_ELEMENTS_INPUT PB0_X0_SLICE_DATA
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
plot_cache="${SLURM_TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID}"
mkdir -p "${plot_cache}"
export MPLCONFIGDIR="${plot_cache}"
python3 py_funcs/plot_pb0_cr33_visualization.py
