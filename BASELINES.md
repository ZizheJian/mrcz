# MRCZ正式方法与baseline口径

本文件固定论文实验中各方法名称的实际含义。`uniform`表示压缩器只看到规则张量的数组下标，不使用GLL物理坐标；`GLL`表示预测或变换显式使用GLL节点位置。所有外部baseline必须直接调用对应压缩器的官方API，自写预测器只作为消融或历史诊断，不称为external baseline。

命令行和bitstream继续使用本文件中的完整名称。CSV及论文图例使用短名称：`MRCZ-3D-GLL`显示为`MRCZ`，`SZ3-3D-GLL-optimized`显示为`FDWZ`，官方API baseline统一省略必然成立的`uniform`，例如`SZ3-4D-uniform`显示为`SZ3-4D`。`FDWZ`表示一种fast distance-weighted interpolation-based lossy compressor，其`Z`作为压缩器名称后缀。

## 论文方法

| 方法 | 定位 | 实际含义 | 当前状态 |
|---|---|---|---|
| `MRCZ-3D-GLL` | 主要方法 | 对每个高阶element做GLL张量基变换、系数与残差量化、32x32分块coefficient-major转置、Huffman和Zstd | 正式方法 |
| `SZ3-3D-GLL-optimized` | 可能的第二个主要方法 | 项目内GLL插值预测器，使用最大有效tree数，并在partition内采样选择维度顺序和linear/cubic | 仅在IPDPS双方法路线成立时进入主要结果 |
| `SZ3-3D-GLL` | 方法消融 | 项目内最简GLL线性插值、单tree、按预测访问顺序保存QB | 不是external baseline |

其他`MRCZ-*`、`SZ3-3D-GLL-*`和`QoZ-3D-GLL`入口均为实现消融、性能原型或诊断工具，不进入正式baseline集合。

## 官方API baseline及tensorization

设一个partition包含`E`个element，单个element的shape为`[N0,N1,N2]`，`N=N0*N1*N2`。表中的维度描述采用逻辑顺序；wrapper会按各官方API规定的最快变化维度传参，但不会重排原始连续数据。

| 方法名模式 | 每个partition的API调用 | 坐标假设 | 额外拼接metadata |
|---|---|---|---|
| `<compressor>-1D-uniform` | 一次一维调用：`[E*N]` | uniform | 无 |
| `<compressor>-2D-uniform` | 一次二维调用：`[E,N]` | uniform | 无 |
| `<compressor>-3D-uniform-fused` | 一次三维调用：`[E*N0,N1,N2]` | uniform | 无 |
| `<compressor>-3D-uniform-per-element` | 每个element一次三维调用：`[N0,N1,N2]` | uniform | 每个partition保存element数量及每条流的`uint64_t`长度 |
| `<compressor>-4D-uniform` | 一次四维调用：`[E,N0,N1,N2]` | uniform | 无 |

各压缩器的实际支持范围如下：

| 压缩器 | 官方API与误差模式 | 已测试布局 | 固定的大数据布局 |
|---|---|---|---|
| SZ3 | SZ3 API，ABS error bound | 1D、2D、fused 3D、per-element 3D、4D | `SZ3-4D-uniform` |
| SZ2 | `SZ_compress_args`，ABS error bound | 1D、2D、fused 3D、per-element 3D、4D | `SZ2-3D-uniform-fused` |
| ZFP | ZFP accuracy mode | 1D、2D、fused 3D、per-element 3D、4D | `ZFP-3D-uniform-fused` |
| SPERR | SPERR C API，pointwise error tolerance | 2D、fused 3D、per-element 3D | `SPERR-2D-uniform` |
| QoZ 2.0 | QoZ q4 API，ABS error bound | 1D、2D、fused 3D、per-element 3D、4D | `QoZ-4D-uniform` |

这些布局完全由小数据实验预先选定，作为pb0和jicf正式实验的固定baseline集合。大数据结果不再用于更换布局，避免根据主要测试数据选择对baseline最有利或最不利的tensorization。

## 公平性与共同封装

- REL输入先使用当前所选field和element范围内的全局`max-min`转换为ABS error bound，再传给后端。
- float输入在`REL<1e-4`时提升到double进行压缩和解压；field bitstream首字节记录实际压缩数据类型。
- 官方SZ3的float路径曾最多超出请求EB约0.031%，因此wrapper传入`0.999*abs_eb`；重新生成的结果均满足用户请求的EB。
- 每个方法使用相同的连续partition边界。除per-element布局外，每个压缩器每个partition只生成一条官方bitstream。
- field文件`<data>_<field>.<method>`包含数据类型、常数标志、partition数量、field最小/最大值、各partition流长度和流本身。BR当前只统计这个field文件，不含共享的`<data>_header.mrcz`与`<data>_element_map.mrcz`。
- 后端计时包含压缩器API内部的预测、量化和熵编码，也包含per-element布局的多次API调用与流拼接；不包含读取原文件、计算field范围、生成共享header/element map、写field文件、写重建文件及误差验证。
- throughput使用上述后端时间计算。论文中若报告端到端时间，必须另行测量并明确标注，不能与backend throughput混用。

## 已知例外

- SZ2会先删除长度为1的维度。对于完整3D element，其4D路径内部合并`E`和`N0`，与fused 3D产生相同bitstream，因此大数据不重复报告两个名称；对于`[N0,N1,1]`的2D element，两者分别退化成`[E,N0,N1]`和`[E*N0,N1]`，仍是不同布局。
- SZ2的官方实现将压缩与解压状态保存在进程级全局变量中，`SZ_decompress_customize_threadsafe`也只是转调普通解压API。完整pb0和jicf上并发调用会产生损坏bitstream或在解压时触发heap corruption，因此正式封装串行调用SZ2 API；其throughput按实际使用的1个CPU报告，不能与支持partition并行的方法作128线程scaling比较。
- MFEM ex1的`5x5x1` per-element SZ2在REL `1e-4`于官方库内部触发allocator断言；该失败保留在原始实验记录中，该布局也已因API及metadata开销淘汰。
- SPERR官方API没有1D和4D入口，不能为了补齐表格而使用自写等价实现。
- QoZ的所有布局都会删除长度为1的退化维度，避免没有物理信息的轴触发官方调优代码错误；例如二维element的per-element布局调用二维API，逻辑4D布局实际调用`[E,N0,N1]`三维API。方法名表示完整三维element时的布局上限。
- QoZ 2.0 q4的自动调优在官方实现中只对二维和三维启用；一维布局仍通过官方QoZ API，但固定使用其插值后端，以避开默认`ALGO_INTERP_LORENZO`在一维上的堆损坏。
- NekRS文件按空间第一维最快存储；读取时将header中的空间shape反转成内部的慢维到快维顺序，因此各种API布局都可直接使用同一段连续数据，不复制或转置field。
- `QoZ-3D-GLL`、`SZ2-3D-uniform`、`SZ2-3D-GLL`和`SZ3-3D-uniform`均为项目内实现，不属于官方API baseline。

## 大数据使用规则

1. 每个官方压缩器保留上表中由小数据选出的一个布局；bitstream完全重复的布局只保留语义最自然的名称。
2. pb0和jicf直接使用这组固定布局生成正式RD与throughput结果，不根据大数据排名重新选择方法。
3. 允许在大数据随机区间上做内部预检查，用于验证程序能运行、误差界限正确、内存与任务时间合理；这些抽样结果不进入论文，也不参与baseline选择。
4. 正式实验固定数据、field、REL EB、element范围、partition数量、CPU数量和计时口径。
5. 将最终方法列表同步到`t1_start_c_code.py`和论文，不再把诊断入口混入正式结果。
