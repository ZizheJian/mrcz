MRCZ experiments
================

# XXX 本目录保存从2026-08-25开始、直接服务于MRCZ和DCC论文的实验脚本及结果。
# XXX SZ3-3D-GLL设计空间实验继续保存在../sz3_3d_gll/，不要混入本目录。
# XXX 大数据只在计算节点运行；计时必须使用程序内部backend时间。

计划保存的结果
==============

- coefficient_eb：全局coeff_eb_scale和模式分组error bound实验。
- lossless_layout：coefficient transpose、符号分布和可能的固定group/tree实验。
- rate_distortion：五个NekRS数据的正式BR-PSNR CSV与作图脚本。
- throughput：不同PSNR和线程数的内部compression/decompression throughput结果。
- ablation：MRCZ最终算法的消融结果。
- regression：float/double、2D/3D、ABS/REL和格式错误测试的简短日志。

文件规则
========

- 脚本使用描述性名称，不继续增加e1、e2式编号。
- Slurm脚本与对应非空.out放在同一子目录；成功任务产生的空.err直接删除。
- CSV只保存可复现的汇总结果，不提交压缩bitstream、解压文件、header或element map临时文件。
- 每组实验的README记录数据路径、field、element范围、error bound、CPU数量、Git commit和结论。
