#!/bin/bash
#SBATCH --job-name=mrcz-e2e-prepare
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

case "${MRCZ_E2E_DATASET:-}" in
    pb0)
        config_name=e2e_pb0_compressed_files
        data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/hyuan/pb0.f00001
        fields=UX,UY,UZ,P,T
        rel_eb=1e-4
        run_cpus=64
        ;;
    jicf0)
        config_name=e2e_jicf0_compressed_files
        data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/viral-flamespeed/jicf0.f00000
        fields=UX,UY,UZ,P,S00
        rel_eb=1e-6
        run_cpus=128
        ;;
    *)
        echo "MRCZ_E2E_DATASET must be pb0 or jicf0" >&2
        exit 2
        ;;
esac

allocated_cpus=${SLURM_CPUS_PER_TASK:-0}
if ((allocated_cpus<run_cpus));then
    echo "This dataset requires at least ${run_cpus} allocated CPUs; got ${allocated_cpus}" >&2
    exit 2
fi

srun --ntasks=1 --cpus-per-task="${run_cpus}" --cpu-bind=cores python3 -u t1_start_c_code.py "${config_name}" --skip-build
python3 py_funcs/create_transfer_benchmark_csv.py "${data_path}" --rel-eb "${rel_eb}" --fields "${fields}"
