#!/bin/bash
#SBATCH --job-name=mrcz-profile
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00
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
method=MRCZ-3D-GLL

cleanup()
{
    rm -f "${data_path}_${field}.${method}" "${data_path}.${method}.bin" "${data_path}_header.mrcz" "${data_path}_element_map.mrcz"
}

trap cleanup EXIT

echo "MRCZ phase profile: data=pb0 field=${field} elements=${elements} REL=${rel_eb} application_cpus=1 allocated_cpus=4"
./install/bin/mrcz_compress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 1 -elements "${elements}" -fields "[${field}]"
./install/bin/mrcz_decompress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 1 -elements "${elements}" -fields "[${field}]"
