#!/bin/bash
#SBATCH --job-name=mrcz-coeff-trees
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=experiments/mrcz_coeff_trees/slurm-%j.out
#SBATCH --error=experiments/mrcz_coeff_trees/slurm-%j.err

set -euo pipefail
cd /home/x-zjian1/mrcz

data_path=/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000
field=S00
methods=(MRCZ-3D-GLL MRCZ-3D-GLL-coeff-3tree MRCZ-3D-GLL-coeff-4tree)
ebs=(1e-4 1e-6 1e-8)
element_ranges=("[0,100000)" "[1500000,1600000)" "[3000000,3100000)")

for elements in "${element_ranges[@]}";do
    for eb in "${ebs[@]}";do
        for method in "${methods[@]}";do
            echo "TEST method=${method} field=${field} REL_EB=${eb} elements=${elements} max_cpus=4 partition_size=25000"
            ./install/bin/mrcz_compress -format nekrs -i "${data_path}" -M REL "${eb}" -method "${method}" -max-cpus 4 -elements "${elements}" -fields "[${field}]"
        done
    done
done

rm -f "${data_path}.MRCZ-3D-GLL.bin" "${data_path}.MRCZ-3D-GLL-coeff-3tree.bin" "${data_path}.MRCZ-3D-GLL-coeff-4tree.bin"
rm -f "${data_path}_${field}.MRCZ-3D-GLL" "${data_path}_${field}.MRCZ-3D-GLL-coeff-3tree" "${data_path}_${field}.MRCZ-3D-GLL-coeff-4tree"
rm -f "${data_path}_header.mrcz" "${data_path}_element_map.mrcz"
