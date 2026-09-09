# ZFP误差分布诊断

Slurm任务`20492785`对pb0五个field在固定约100 dB设置下的ZFP重建误差进行了全量统计，共检查4,984,714,240个值。所有field均没有误差达到`0.5*ABS_EB`，因此也没有误差达到或超过ABS EB；最大的归一化绝对误差来自T，为`0.21956094*ABS_EB`。各field的中位误差不超过`0.0161*ABS_EB`，99%误差不超过`0.0700*ABS_EB`，说明该设置下ZFP整体显著低于名义误差上限，而不是仅有少量大误差触及上限。

逐field的实际PSNR、最大误差、归一化误差分位数和阈值计数保存在`zfp_pb0_100db_summary.csv`。统计由解压命令的`-analyze-error-distribution true`开关执行，位于backend解压计时段之后，因此不影响正式decompression time；默认关闭时不执行直方图统计。
