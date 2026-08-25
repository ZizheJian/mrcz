#!/bin/bash
#SBATCH --job-name=mrcz-simd-4cpu
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
export OMP_NUM_THREADS=4
export MRCZ_PROFILE_PHASES=1

data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000
field=S00
method=MRCZ-3D-GLL
ebs=(1e-4 1e-6 1e-8)
element_ranges=("[0,100000)" "[1500000,1600000)" "[3000000,3100000)")

cleanup()
{
    rm -f "${data_path}_${field}.${method}" "${data_path}.${method}.bin" "${data_path}_header.mrcz" "${data_path}_element_map.mrcz"
}

trap cleanup EXIT

for elements in "${element_ranges[@]}";do
    for rel_eb in "${ebs[@]}";do
        echo "MRCZ SIMD multicore profile: field=${field} elements=${elements} REL=${rel_eb} application_cpus=4 partition_size=25000"
        ./install/bin/mrcz_compress -format nekrs -i "${data_path}" -M REL "${rel_eb}" -method "${method}" -max-cpus 4 -elements "${elements}" -fields "[${field}]"
        cleanup
    done
done
