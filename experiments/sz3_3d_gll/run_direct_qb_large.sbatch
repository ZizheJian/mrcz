#!/bin/bash
#SBATCH --job-name=sz3-direct-qb
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00
#SBATCH --output=experiments/sz3_3d_gll/slurm-%x-%j.out
#SBATCH --error=experiments/sz3_3d_gll/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1

data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001
field=T
elements='[1000000,1100000)'
rel_eb=1e-4

cleanup_method()
{
    local method=$1
    rm -f "${data_path}_${field}.${method}" "${data_path}.${method}.bin"
}

cleanup_all()
{
    cleanup_method SZ3-3D-GLL
    cleanup_method SZ3-3D-GLL-optimized
    rm -f "${data_path}_header.mrcz" "${data_path}_element_map.mrcz"
}

trap cleanup_all EXIT

for method in SZ3-3D-GLL SZ3-3D-GLL-optimized
do
    echo "===== ${method} compression ====="
    ./install/bin/mrcz_compress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 1 -elements "${elements}" -fields "[${field}]"
    echo "===== ${method} decompression ====="
    ./install/bin/mrcz_decompress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 1 -elements "${elements}" -fields "[${field}]"
    cleanup_method "${method}"
done
