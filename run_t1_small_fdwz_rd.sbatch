#!/bin/bash
#SBATCH --job-name=mrcz-small-fdwz-rd
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --hint=nomultithread
#SBATCH --time=00:10:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MRCZ_T1_MAX_CPUS=1
export MRCZ_T1_BUILD_BEFORE_TEST=false
export MRCZ_T1_METHODS=FDWZ
export MRCZ_T1_EBS=1e-2,1e-3,1e-4,1e-5,1e-6
export MRCZ_T1_FDWZ_SAMPLE_INTERVAL=128
for dataset in mixlay0 test0 tc_channel0 mfem_ex0 mfem_ex1 nektar_cylinder
do
    export MRCZ_T1_DATASET="${dataset}"
    srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores python3 -u t1_start_c_code.py
done
python3 py_funcs/merge_method_csv.py third_party/pySEMTools/examples/data results/RD/mixlay0 --method FDWZ --pattern 'mixlay0.f00001_*.csv'
python3 py_funcs/merge_method_csv.py third_party/fld_pyutils/demos/data results/RD/test0 --method FDWZ --pattern 'test0.f00001_*.csv'
python3 py_funcs/merge_method_csv.py third_party/pySEMTools/examples/data results/RD/tc_channel0 --method FDWZ --pattern 'tc_channel0.f00001_*.csv'
python3 py_funcs/merge_method_csv.py build/mfem_build/examples/ex0examples results/RD/mfem_ex0 --method FDWZ --pattern 'sol.gf_*.csv'
python3 py_funcs/merge_method_csv.py build/mfem_build/examples/ex1examples results/RD/mfem_ex1 --method FDWZ --pattern 'sol.gf_*.csv'
python3 py_funcs/merge_method_csv.py third_party/flow-stability-cylinder/tutorial/base results/RD/nektar_cylinder --method FDWZ --pattern 'Cylinder-Base.vtu_*.csv'
python3 py_funcs/plot_small_rd.py
python3 py_funcs/summarize_rd.py
