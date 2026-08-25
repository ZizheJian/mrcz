#!/bin/bash
#SBATCH --job-name=mrcz-batched
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:05:00
#SBATCH --output=experiments/mrcz_simd/slurm-%x-%j.out
#SBATCH --error=experiments/mrcz_simd/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1
export MRCZ_PROFILE_PHASES=1

data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001
field=T
elements='[1000000,1100000)'
rel_eb=1e-4
methods=(MRCZ-3D-GLL MRCZ-3D-GLL-element-batched)

cleanup()
{
    for method in "${methods[@]}";do
        rm -f "${data_path}_${field}.${method}" "${data_path}.${method}.bin"
    done
    rm -f "${data_path}_header.mrcz" "${data_path}_element_map.mrcz"
}

trap cleanup EXIT

for repeat in 1 2 3;do
    for method in "${methods[@]}";do
        echo "MRCZ element-batched large comparison: repeat=${repeat} method=${method} field=${field} elements=${elements} REL=${rel_eb} application_cpus=1"
        ./install/bin/mrcz_compress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 1 -elements "${elements}" -fields "[${field}]"
        ./install/bin/mrcz_decompress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 1 -elements "${elements}" -fields "[${field}]"
    done
done

cmp "${data_path}_${field}.MRCZ-3D-GLL" "${data_path}_${field}.MRCZ-3D-GLL-element-batched"
cmp "${data_path}.MRCZ-3D-GLL.bin" "${data_path}.MRCZ-3D-GLL-element-batched.bin"
echo "Large-sample bitstreams are byte-identical."
