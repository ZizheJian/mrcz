MRCZ partition-local QB实验
===========================

# XXX 日期：2026-08-25。
# XXX 对照版本：独立partition_num已经实现，但coefficient QB和residual QB仍各自覆盖整个field。
# XXX 新版本：每个OpenMP线程各分配一次足以容纳其最大partition的两份未初始化int64_t数组，随后循环复用。

正确性
======

- test0/T、REL 1e-3、8个partition：正式MRCZ的新旧bitstream SHA-256均为
  5024e30863982b73b2aca211acbb87a666bd2228f67a6c2439615390ee66ebb8。
- test0/T、REL 1e-5触发float到double提升，4 CPU和8个partition下新旧bitstream逐字节相同。
- MRCZ-3D-GLL-element-batched和MRCZ-3D-GLL-element-batched-original-transpose使用7个非整batch partition时，
  新旧bitstream也逐字节相同。
- 正式MRCZ解压的BR、PSNR、MSE和最大误差不变。

单核内存与速度
==============

输入为pb0/T前50000个element，REL 1e-3，固定8个partition和1 CPU。输入通过/tmp中的符号链接读取，
内部Compression time/Decompression time用于速度比较，/usr/bin/time -v只用于读取Maximum resident set size。
新旧版本交错运行5次，以下为中位数：

                         全field QB       partition-local QB      变化
Compression time         0.600063 s       0.578747 s              -3.55%
Compression peak RSS     594920 KB        244808 KB               -58.85%
Decompression time       0.369172 s       0.340020 s              -7.90%
Decompression peak RSS   540568 KB        240596 KB               -55.49%

多核初步结果
============

同一输入使用4 CPU和8个partition时，压缩峰值RSS中位数由818568降至618620 KB，解压由620976降至
421044 KB。登录节点上的压缩时间噪声较大，新版本的5次中位数约慢5.8%，而解压约快16.3%；这两个速度变化方向相反，
不能作为正式结论。后续应在独占计算节点、关闭阶段profiling并交错运行至少7次。

理论内存
========

jicf单个10x10x10 field共有3240264000个值。两份全field int64_t QB占51844224000 bytes，即约51.84 GB
或48.28 GiB。固定256个partition时，partition-local两份QB在1线程下约为0.203 GB，在128线程下约为25.92 GB。
Huffman临时缓冲区仍然与每个活跃partition的两份QB近似等大，因此1线程float压缩的估算峰值可由约65.0 GB降至
约13.4 GB；严格EB自动提升double时可由约78.0 GB降至约26.3 GB。解压后的完整field和用于误差检查的原始field
仍可能成为下一阶段的内存瓶颈。
