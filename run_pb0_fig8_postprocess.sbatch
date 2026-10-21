#!/bin/bash
#SBATCH --job-name=fig8-plot
#SBATCH --account=cis240669
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=results/RD/slurm-%x-%j.out
#SBATCH --error=results/RD/slurm-%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /home/x-zjian1/miniconda3/etc/profile.d/conda.sh
conda activate "${SLURM_SUBMIT_DIR}/mrcz_env"
export MPLCONFIGDIR="/tmp/mrcz-matplotlib-${SLURM_JOB_ID}"
python3 py_funcs/plot_fig8_transform_comparison.py
cd 论文
/apps/spack/anvil/apps/texlive/20200406-gcc-11.2.0-eeavxnm/bin/x86_64-linux/pdflatex -interaction=nonstopmode -halt-on-error main.tex
/apps/spack/anvil/apps/texlive/20200406-gcc-11.2.0-eeavxnm/bin/x86_64-linux/pdflatex -interaction=nonstopmode -halt-on-error main.tex
