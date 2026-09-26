# 固定REL EB端到端实验

pb0使用`REL EB=1e-4`和64 CPU，jicf0使用`REL EB=1e-6`和128 CPU；两者都是正式五点RD扫描在对数REL EB轴上的中点。压缩方法将五个逐field bitstream作为一个Globus任务传输，Direct则传输一个包含五个目标field的原始数据文件。每种配置均从ACCESS Anvil传输到LCRC Improv DTN并独立重复三次。固定PSNR的旧实验文件和传输CSV仅作为历史结果保留，不参与当前汇总。

- `pb0_transfer_rel_eb.csv`和`jicf0_transfer_rel_eb.csv`：数据目录旁原始测速CSV的版本化快照。实际速度分别填写到`hyuan/pb0.f00001_transfer_rel_eb_0.0001.csv`和`viral-flamespeed/jicf0.f00000_transfer_rel_eb_1e-06.csv`，再同步到这里用于论文汇总。七个压缩方法的每次任务包含五个逐field bitstream，Direct的每次任务包含一个原始数据文件；每种配置均有三行trial。
- `end_to_end_summary.csv`：五个field的RD backend压缩/解压时间、三次传输时间、压缩文件大小、CR和端到端排名。
- `end_to_end_total.svg`：压缩、三次传输的中位数和解压时间组成的堆叠柱状图；图中不绘制误差棒。
- `slurm-mrcz-e2e-*-rel-*.out`：固定REL EB压缩文件生成任务`20493296`和`20493297`的标准输出；两个任务的标准错误均为空，因此不保留空`.err`文件。
- `run_end_to_end_rel_eb.sbatch`：只执行压缩、保留默认位置输出且不改写RD CSV的Slurm入口；不申请exclusive或显式内存。
- `direct_transfer_manifest.csv`：两个Direct文件的field顺序、路径、字节数和SHA-256；由`run_prepare_direct_transfer.sbatch`生成。

Direct文件位于`/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/direct_transfer/{pb0,jicf0}`。每个文件是五个目标field在原始NekRS文件中的连续数据段，不包含坐标、未评估field或文件头。每次重复应将目标端对应文件删除或换用新的目标目录，以免Globus跳过已经存在且校验一致的文件。

运行`python3 py_funcs/summarize_end_to_end.py`重新生成汇总CSV，再运行`python3 py_funcs/plot_end_to_end.py`生成SVG。汇总脚本使用`Bytes Transferred/effective_bytes_per_second`计算每次传输时间，并以per-CPU MiB/s乘申请CPU数还原并行backend时间。Direct的速度留空时，汇总会保留一个总时间为0的占位项；填入三次速度后会自动替换为实际传输时间。SZ2不进入端到端汇总，因为其官方API使用进程级可变全局状态，无法安全地参与多核工作流比较。
