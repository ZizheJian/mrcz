# 1. 简介

本项目研究高阶有限元与谱元模拟数据的误差有界有损压缩，主要面向NekRS/Nek5000等程序生成的三维物理场，并使用MFEM和Nektar++样例验证数据布局的通用性。此类数据由大量高阶element组成：element之间按非规则网格及拓扑关系组织，而每个element内部的采样点通常位于Gauss--Lobatto--Legendre（GLL）节点上，形成尺寸较小但结构固定的张量积网格。物理场在单个element内部通常较平滑，但跨element的存储顺序、尺度和局部变化可能不连续，因此把数据简单展平或视为普通规则网格会损失可利用的结构。

项目的核心目标是在严格满足逐点误差界限的前提下，同时改善压缩率、压缩与解压吞吐量以及多核扩展性。主要方法MRCZ对每个element进行GLL张量基变换，分别编码低频系数与点值残差，再通过coefficient-major布局、Huffman和Zstd压缩；候选方法FDWZ则直接利用GLL节点距离进行自适应插值预测。两者将与使用官方API和预先固定tensorization方式的通用压缩器进行公平比较，并以rate--distortion、throughput、scaling和端到端传输时间共同决定最终论文定位。

项目主入口为`t1_start_c_code.py`，其中定义正式数据集、field、误差界限、方法列表以及结果文件位置。正式CSV保存在对应数据旁边，方法与baseline口径见`BASELINES.md`；本文件只记录投稿IPDPS或DCC之前仍需完成的工作。

# 2. FDWZ核心与优化记录

FDWZ（Fast Distance-Weighted Interpolation-Based Lossy Compressor）借鉴SZ3的多层插值预测框架，但预测范围限制在单个高阶element内，并使用真实GLL节点坐标而不是均匀索引距离。它不是对官方SZ3的简单参数调整，而是项目内为高阶有限元/谱元数据实现的独立预测器。

## 2.1. 核心流程

1. 根据element shape和各维GLL坐标预先构造可由编解码两端确定性复现的预测层、访问顺序、reference位置及插值权重。每个有效维度的两个端点构成初始边界，递归选择最接近物理中点的未编码GLL节点。
2. 每个二维/三维element的4/8个角点直接进行误差有界量化；其余节点按照预测层和维度顺序逐步处理。linear模式使用两个已重建reference按真实GLL距离进行一次插值；cubic模式则根据当前可用的已重建reference数量，自适应使用三次、二次或一次Lagrange插值。
3. 每次预测后立即量化prediction residual，并用量化结果覆盖目标值。后续预测只引用该重建值，因此压缩端与解压端具有完全相同的误差传播路径，并保证逐点误差界限。
4. FDWZ在每个partition内进行确定性的等间隔采样，通常采样128--1024个element，并取各等分区间中点附近的element；当partition不足128个element时使用全部element。随后先从六种维度顺序中选择sample真实压缩大小最小的linear方案，再在相同顺序下比较linear与cubic。选中的方案编号写入该partition的bitstream。
5. 量化整数按照角点及预测层分组。正式FDWZ固定使用当前element shape允许的最大Huffman tree数，最多五棵；每个partition独立建立Huffman tree，将全部Huffman metadata和payload拼接后再压缩为一个Zstd frame。
6. partition之间彼此独立，通过OpenMP并行执行采样、预测、量化、Huffman和Zstd；解压端按相反顺序恢复Huffman符号，并在每个element内依次重建节点。

## 2.2. 角点与初始seed阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| 4/8个角点预测为0并直接量化 | **采纳** | 最简单且在五个NekRS数据上最稳定，作为正式seed方案。 |
| 使用前一编号element预测角点 | **放弃** | 依赖编号相邻element具有物理相关性，完整bitstream没有稳定收益。 |
| element内角点Lorenzo、角点均值加deviation及其组合 | **放弃** | 局部角点流偶有改善，但额外符号或后续误差传播抵消了完整BR收益。 |
| 将角点流从element-contiguous转成local-symbol-contiguous | **放弃** | test0、mixlay0和jicf0退化，仅tc_channel0和pb0小幅改善。 |
| 先建立中心、最后向外外插角点 | **放弃** | 角点/中心点预测误差诊断已排除该方向。 |

## 2.3. element内预测与访问顺序阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| 按GLL物理中点递归并使用真实GLL距离计算权重 | **采纳** | 比索引中点更符合非均匀GLL节点分布，是正式预测树的基础。 |
| 三维数据中的六种维度顺序 | **采纳为partition级选择** | 不存在跨数据和partition普遍最优的固定顺序，因此由确定性的等间隔sample选择。（SZ3中已有） |
| linear与cubic | **采纳为partition级选择** | 高阶插值在mixlay0、tc_channel0、pb0和部分jicf0上显著降低BR，但会破坏test0的特殊可预测结构，不能全局固定启用。（SZ3中已有） |
| surface-first：按边、面、体及中心向外顺序预测 | **放弃** | 在代表性数据上没有形成rate--distortion优势。 |
| normal-boundary：先建立中心块，再沿法向外插边界 | **放弃** | 五个数据上均被正式物理中点预测树Pareto支配。 |
| 官方SZ3式2幂stride、单原点seed的index-order | **放弃** | 对小element反复支付粗层级远距离预测成本，三个小数据均无优势。 |
| 不预测物理中点，而是预测索引中点 | **放弃** | pb0和jicf0的BR分别增加约8.4%和4.7%，说明物理距离比索引邻近更重要。 |
| element编号方向预测 | **放弃** | 在五个数据上严格优于所有element内方向的点不足约1.1%，加入空间oracle后的总体改善很小。 |
| 真实共享面邻居与坐标匹配预测 | **放弃当前实现方向** | 早期跨公共面实验没有优于element内预测，同时会引入邻居搜索和拓扑metadata成本。 |
| QoZ式GLL预测路径 | **仅保留独立对照** | test0/UZ上的BR明显差于FDWZ，不作为正式优化继续。 |

## 2.4. 量化与residual建模阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| 对之后仍会被引用的reference使用更严格EB | **放弃** | 同PSNR下仅mixlay0约有2%收益，其余数据持平或退化，不足以承担选择复杂度。 |
| 每个局部位置学习prediction-residual bias | **放弃** | 单bin峰值和带置信度的稳健宽窗口版本都不能覆盖采样与bias metadata开销，且降低吞吐量。 |
| 根据前面element在线更新预测统计 | **放弃** | 会引入顺序依赖和同步成本，没有进入正式bitstream。 |

## 2.5. QB布局、Huffman与Zstd阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| 非角点QB按实际预测访问顺序保存 | **采纳** | 大数据随机段上稳定降低约5%--7% BR，压缩吞吐代价约1%，解压局部性也更好。（SZ3中已有） |
| 非角点QB转置一下，使得element编号维度最快 | **放弃** | pb0和jicf0平均BR分别增加约4.2%和3.5%，跨步访问还显著破坏缓存局部性。 |
| 1--5棵Huffman tree分离角点和不同预测层 | **采纳最大有效tree数** | 大partition上多tree可利用各预测层不同的符号分布；正式FDWZ按shape固定最多五棵，避免额外tree选择逻辑。 |
| 从sample直接选择tree数 | **放弃** | 小sample无法正确摊销tree metadata，tree数与完整partition的最优选择经常不一致。 |
| 按符号数阈值或half/full sample外推tree成本 | **放弃** | regret仍不够稳定；最终固定最大有效tree比自适应规则更简单稳健。 |
| 每个partition的多个Huffman流合并到一个Zstd frame | **采纳** | 保留不同符号分布的Huffman收益，同时让Zstd继续利用整个partition的字节级冗余。 |

## 2.6. 配置采样、并行与实现加速阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| partition内确定性的等间隔采样 | **采纳** | 对维度顺序和linear/cubic的holdout选择较稳定，最多采样1024个element，不依赖随机数或随机种子。 |
| 只采样partition连续前缀 | **放弃** | 对非平稳partition代表性不足，曾出现超过100%的holdout regret。 |
| 预计算prediction stencil和GLL权重 | **采纳** | 避免在每个element重复构造预测关系，并保证编解码路径一致。 |
| partition级OpenMP并行 | **采纳** | partition天然独立，无需在element预测树内部增加同步。 |
| 跨element AVX2并行预测与量化 | **放弃正式合并** | 约256 KB tile解决了整partition position-major访问的缓存退化，但float没有收益，double端到端仅约快9%，不足以承担独立舍入语义和复杂实现。 |

上述实验的核心结论记录在本节。正式FDWZ对应代码方法名为`SZ3-3D-GLL-optimized`，论文和CSV图例使用`FDWZ`。

说明：2026-08-29以前的详细算法与实验记录见experiments/TODO_history_2026-08-24.txt和
experiments/TODO_history_2026-08-29.txt。本文件只保留投稿前仍需完成的工作。

# 3. 当前状态

- 正式MRCZ算法已冻结：逐element DChT、系数与残差量化、32x32分块coefficient-major
  转置、Huffman和Zstd，partition之间OpenMP并行。
- FDWZ已冻结：最大有效Huffman tree数，并在每个partition采样选择维度顺序和
  linear/cubic预测。
- 官方API baseline及tensorization已经冻结，定义见BASELINES.md。
- pb0和jicf0的完整rate-distortion与128 CPU throughput CSV已经生成，保存在数据旁边。
- pb0 100 dB与jicf0 150 dB的fixed-quality scaling实验已经完成，结果见
  experiments/mrcz_scaling/。

# 4. 投稿前主线

[ ] 1. 冻结实验版本与结果。
       提交当前正式代码、实验脚本和结果说明；记录生成正式CSV所用的commit、Slurm job、
       编译器、CPU、partition数量和计时口径。CSV继续保存在对应数据旁边。

[ ] 2. 生成论文用rate-distortion汇总。
       从pb0和jicf0各field CSV生成BR--PSNR图，并在固定PSNR处汇总MRCZ、FDWZ和五个
       官方baseline的BR、CR、最大误差与异常点。旧稿中的数字不得直接复用。

[ ] 3. 整理速度与扩展性结果。
       绘制single-core throughput、fixed-quality strong/weak scaling、speedup和parallel
       efficiency；明确64--128 CPU后的内存带宽与NUMA瓶颈，不声称近线性扩展。

[ ] 4. 完成新版端到端实验。
       重新测量direct、MRCZ和FDWZ的压缩--传输--解压总时间，并使用
       T_total=T_compress+compressed_size/bandwidth+T_decompress扫描存储或网络带宽，
       给出三种方案各自占优的区间。旧稿Globus结果只作历史参考。

[ ] 5. 决定投稿方向。
       若MRCZ与FDWZ在rate--throughput--end-to-end空间形成稳定互补Pareto点，按IPDPS
       Algorithms组织；否则按DCC组织，只把MRCZ作为主要算法。

[ ] 6. 重写论文。
       使用high-order finite/spectral-element field data描述研究对象；替换旧baseline、
       实验数字、方法名称和贡献列表。准备数据布局与流程图、BR--PSNR图、scaling图、
       消融或bitstream组成图，以及数据集汇总表。

# 5. 投稿后或artifact阶段

[ ] 7. 审计并版本化bitstream metadata，增加未知版本、截断文件、错误method和错误数据类型测试。

[ ] 8. 让t1在超过error bound、重建数量错误或编解码失败时自动返回非零状态，并处理低于机器精度的EB。

[ ] 9. 清理正式路由中剩余的实验方法和历史注释；完成安装后header-only target及CLI复现检查。

# 6. 不进入当前论文

- 坐标X/Y/Z专用有损压缩及几何/拓扑安全定义。
- NekRS尾metadata的语义解析；当前继续作为opaque blob无损保存。
- coefficient delta、更多SZ3预测方案、跨element AVX2和GPU实现。
