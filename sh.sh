#!/bin/bash
#SBATCH -J myjob                # 作业名
#SBATCH -N 1                    # 1 个节点
#SBATCH -n 1                    # 1 个任务
#SBATCH -c 1                    # 每个任务 1 核心
#SBATCH --time=1:00:00          # 时间（可改）
#SBATCH --output=log.out        # 输出日志
#SBATCH --error=log.err         # 错误日志

source ~/.bashrc
conda activate unsz_env/
python3 t1_convert_nek5000_to_binary.py