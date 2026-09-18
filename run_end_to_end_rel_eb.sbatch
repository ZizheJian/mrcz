#!/bin/bash
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --hint=nomultithread
#SBATCH --time=00:30:00
#SBATCH --output=results/end_to_end/slurm-%x-%j.out
#SBATCH --error=results/end_to_end/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export MRCZ_T1_BUILD_BEFORE_TEST=false
export MRCZ_T1_METHODS=MRCZ,FDWZ,SZ3-4D,SZ2-3D-fused,QoZ-4D,ZFP-3D-fused,SPERR-2D
export MRCZ_T1_CLEANUP_AFTER_EACH_TEST=false
export MRCZ_T1_RUN_DECOMPRESSION=false
export MRCZ_T1_WRITE_TO_CSV=false
export MRCZ_T1_FDWZ_SAMPLE_INTERVAL=128

case "${MRCZ_E2E_DATASET}" in
    pb0)
        data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001
        export MRCZ_T1_DATASET=pb0
        export MRCZ_T1_MAX_CPUS=64
        export MRCZ_T1_FIELDS=UX,UY,UZ,P,T
        export MRCZ_T1_EBS=1e-4
        ;;
    jicf0)
        data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000
        export MRCZ_T1_DATASET=jicf0
        export MRCZ_T1_MAX_CPUS=128
        export MRCZ_T1_FIELDS=UX,UY,UZ,P,S00
        export MRCZ_T1_EBS=1e-6
        ;;
    *)
        echo "MRCZ_E2E_DATASET must be pb0 or jicf0" >&2
        exit 2
        ;;
esac

srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores python3 -u t1_start_c_code.py
python3 py_funcs/create_transfer_csv.py "${data_path}" --rel-eb "${MRCZ_T1_EBS}" --fields "${MRCZ_T1_FIELDS}"
