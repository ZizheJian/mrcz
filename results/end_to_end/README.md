# 固定REL EB端到端实验

pb0使用`REL EB=1e-4`和64 CPU，jicf0使用`REL EB=1e-6`和128 CPU；两者都是正式五点RD扫描在对数REL EB轴上的中点。每种方法将五个field文件作为一个Globus任务从ACCESS Anvil传输到LCRC Improv DTN，并独立重复三次。固定PSNR的旧实验文件和传输CSV仅作为历史结果保留，不参与当前汇总。

- `pb0_transfer_rel_eb.csv`和`jicf0_transfer_rel_eb.csv`：数据目录旁原始传输CSV的版本化快照，速度来自Globus debug data的`effective_bytes_per_second`。
- `end_to_end_summary.csv`：五个field的RD backend压缩/解压时间、三次传输时间、压缩文件大小、CR和端到端排名。
- `end_to_end_total.svg`：端到端时间中位数；误差棒仅表示三次传输带来的最小--最大范围。
- `slurm-mrcz-e2e-*-rel-*.out`：固定REL EB压缩文件生成任务`20493296`和`20493297`的标准输出；两个任务的标准错误均为空，因此不保留空`.err`文件。
- `run_end_to_end_rel_eb.sbatch`：只执行压缩、保留默认位置输出且不改写RD CSV的Slurm入口；不申请exclusive或显式内存。

运行`python3 py_funcs/summarize_end_to_end.py`重新生成汇总CSV，再运行`python3 py_funcs/plot_end_to_end.py`生成SVG。汇总脚本使用`Bytes Transferred/effective_bytes_per_second`计算每次传输时间；并行方法以per-CPU MiB/s乘申请CPU数还原backend时间，SZ2因官方API全局状态而固定按1个活动CPU计算。
