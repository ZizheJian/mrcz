# MRCZ

## 脚本列表

| 脚本 | 功能 |
|---|---|
| `t1_start_c_code.py` | 构建并运行C++压缩与解压程序，执行主要rate-distortion和throughput实验，并将各field的结果写入CSV。数据集、方法、field和error bound可通过文件配置或`MRCZ_T1_*`环境变量选择。 |
| `t2_find_three_irregular_elements.py` | 在NekRS数据中寻找三个两两共享面的非平面element，并生成element、内部节点层和reference element的交互式三维图。 |
| `t3_test_interpolation_prediction.py` | 抽取NekRS element、寻找共享面邻居，比较跨公共面预测与element内部预测的误差，并生成相对误差ECDF。 |
| `t4_plot_data_coeff_err_distribution.py` | 比较原始field、变换系数和误差的概率分布，输出三联直方图；运行前需要准备对应的`_coeff`和`_err`二进制文件。 |
| `t5_start_weak_scaling.py` | 在固定目标PSNR下执行单节点weak scaling实验，结果写入`results/weak_scaling/`中的raw、summary和quality CSV。 |
| `generate_mfem_and_nektar_data.py` | 构建MFEM ex0/ex1和Nektar++ FieldConvert，并生成项目使用的两个MFEM样例及Nektar++ cylinder VTU数据。 |

正式weak scaling可分别使用`run_t5_weak_pb0.sbatch`和`run_t5_weak_jicf.sbatch`提交到Slurm。

`py_funcs/`中的辅助功能：

| 文件 | 功能 |
|---|---|
| `read_nekrs_header.py` | 解析NekRS field文件头，返回字节序、数据精度、element shape、element数量和field信息。 |
| `interpolate_br_at_psnr.py` | 整理BR--PSNR曲线，并线性插值得到指定PSNR下的BR。 |
| `plot_small_rd.py` | 根据18个小数据field的正式CSV生成RD汇总图。 |
| `summarize_rd.py` | 对每个field在七种正式方法共同PSNR区间的近中点处汇总CR增长，并生成逐field CSV。 |
| `plot_large_rd.py` | 为pb0或jicf0单独生成五个field的正式RD图。 |
| `plot_large_throughput.py` | 根据pb0和jicf0的正式CSV生成以REL EB或PSNR为横坐标的压缩与解压throughput图。 |
| `plot_weak_scaling.py` | 根据weak scaling summary CSV生成`results/weak_scaling/weak_scaling_efficiency.svg`。 |
| `plot_weak_scaling_dataset.py` | 为pb0或jicf0单独生成weak scaling efficiency图。 |

## 安装

需要Python 3、CMake、支持C++20和OpenMP的C++编译器。绘图与分析脚本还需要NumPy、SciPy、Matplotlib和Plotly。

克隆项目并初始化submodule：

```bash
git clone --recurse-submodules git@github.com:ZizheJian/mrcz.git
cd mrcz
```

如果项目已经克隆：

```bash
git submodule update --init --recursive
```

创建Python环境：

```bash
python3 -m venv mrcz_env
source mrcz_env/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib plotly
```

推荐直接运行`t1_start_c_code.py`。它会依次构建并安装SPERR、MFEM和MRCZ，然后执行当前配置的实验：

```bash
python3 t1_start_c_code.py
```

FDWZ默认在每个partition均匀采样六种linear维度顺序，再以最佳linear的维度顺序执行正式cubic压缩。默认采样间隔为128个element且每个partition至少取32个均匀样本；命令行可通过`-fdwz-sample-interval N`和`-fdwz-min-sample-num N`修改，`t1_start_c_code.py`对应使用`MRCZ_T1_FDWZ_SAMPLE_INTERVAL`和`MRCZ_T1_FDWZ_MIN_SAMPLE_NUM`。`-fdwz-selector`还可切换`linear-only`或`cubic-exhaustive`。

安装后的可执行文件位于`install/bin/mrcz_compress`和`install/bin/mrcz_decompress`。已有安装或多个实验并行时，可设置`MRCZ_T1_BUILD_BEFORE_TEST=false`跳过重复构建。

如果只想安装而不运行实验，或者需要单独控制CMake配置，也可以手动构建：

```bash
cmake -S third_party/sperr -B build/sperr_build -DCMAKE_INSTALL_PREFIX="$PWD/install/sperr_install" -DBUILD_SHARED_LIBS=OFF
cmake --build build/sperr_build --parallel 8
cmake --install build/sperr_build

cmake -S third_party/mfem -B build/mfem_build -DCMAKE_INSTALL_PREFIX="$PWD/install/mfem_install"
cmake --build build/mfem_build --parallel 8
cmake --install build/mfem_build

cmake -S c_codes -B build -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel 8
cmake --install build
```

安装后，外部C++项目可以通过CMake导入header-only MRCZ核心：

```cmake
find_package(MRCZ CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE mrcz::mrcz)
```

程序包含公共header：

```cpp
#include <mrcz_core/mrcz.hpp>
```

如果MRCZ未安装到系统默认路径，配置外部项目时通过`CMAKE_PREFIX_PATH`指定MRCZ安装目录；MRCZ会自动传递C++20、OpenMP、Zstd、tinyxml2和MFEM依赖。

submodule中提供了三个可直接测试的NekRS小数据：

- `third_party/pySEMTools/examples/data/mixlay0.f00001`
- `third_party/pySEMTools/examples/data/tc_channel0.f00001`
- `third_party/fld_pyutils/demos/data/test0.f00001`

MFEM和Nektar++验证数据可通过以下命令重新生成：

```bash
python3 generate_mfem_and_nektar_data.py
```

脚本生成`build/mfem_build/examples/ex0examples/`、`build/mfem_build/examples/ex1examples/`和`third_party/flow-stability-cylinder/tutorial/base/Cylinder-Base.vtu`。

# 1. 简介

本项目研究高阶有限元与谱元模拟数据的误差有界有损压缩，主要面向NekRS/Nek5000等程序生成的三维物理场，并使用MFEM和Nektar++样例验证数据布局的通用性。此类数据由大量高阶element组成：element之间按非规则网格及拓扑关系组织，而每个element内部的采样点通常位于Gauss--Lobatto--Legendre（GLL）节点上，形成尺寸较小但结构固定的张量积网格。物理场在单个element内部通常较平滑，但跨element的存储顺序、尺度和局部变化可能不连续，因此把数据简单展平或视为普通规则网格会损失可利用的结构。

项目的核心目标是在严格满足逐点误差界限的前提下，同时改善压缩率、压缩与解压吞吐量以及多核扩展性。主要方法MRCZ对每个element进行GLL张量基变换，分别编码变换系数与点值残差，再通过coefficient-major布局、Huffman和Zstd压缩；候选方法FDWZ则直接利用GLL节点距离进行自适应插值预测。两者将与使用官方API和预先固定tensorization方式的通用压缩器进行公平比较，并以rate--distortion、throughput、scaling和端到端传输时间共同决定最终论文定位。

所有方法共用field级外围流程，包括读取数据、计算值域与ABS/REL error bound、选择计算精度以及写入公共bitstream信息。若一个field的最大值和最小值可以由其中点在`abs_eb`内同时重建，公共层会只保存constant flag和该中点，并跳过具体backend；因此constant compression不是任何单个算法的组成部分，可通过`param.if_constant_compression`控制。

项目主入口为`t1_start_c_code.py`，其中定义正式数据集、field、误差界限、方法列表以及结果文件位置。正式CSV保存在对应数据旁边，方法与baseline布局见第5章；本文档后半部分记录投稿IPDPS或DCC之前仍需完成的工作。

# 2. FDWZ核心与优化记录

FDWZ（Fast Distance-Weighted Interpolation-Based Lossy Compressor）借鉴SZ3的多层插值预测框架，但预测范围限制在单个高阶element内，并使用真实GLL节点坐标而不是均匀索引距离。它不是对官方SZ3的简单参数调整，而是项目内为高阶有限元/谱元数据实现的独立预测器。

## 2.1. 核心流程

1. 根据element shape和各维GLL坐标预先构造可由编解码两端确定性复现的预测层、访问顺序、reference位置及插值权重。每个有效维度的两个端点构成初始边界，递归选择最接近物理中点的未编码GLL节点。
2. 每个二维/三维element的4/8个角点直接进行误差有界量化；其余节点按照预测层和维度顺序逐步处理。linear模式使用两个已重建reference按真实GLL距离进行一次插值；cubic模式则根据当前可用的已重建reference数量，自适应使用三次、二次或一次Lagrange插值。
3. 每次预测后立即量化prediction residual，并用量化结果覆盖目标值。后续预测只引用该重建值，因此压缩端与解压端具有完全相同的误差传播路径，并保证逐点误差界限。
4. FDWZ在每个partition内进行确定性的均匀采样，默认采样间隔为128且至少取32个样本。sample真实压缩六种linear维度顺序，再以最佳linear的方向执行正式cubic压缩；选中的完整方案编号写入该partition的bitstream。
5. 量化整数按照角点及预测层分组。正式FDWZ固定使用当前element shape允许的最大Huffman tree数，最多五棵；每个partition独立建立Huffman tree，将全部Huffman metadata和payload拼接后再压缩为一个Zstd frame。QB固定使用32位存储并保留`INT32_MIN`作为unpredictable标记；超出范围或无法严格满足误差界的点写入标记，并在全部Huffman流后按出现顺序追加原始浮点值。残差QB优先由relative EB、field range和全部stencil权重在field级证明安全，只有无法证明时才逐点检查；角点因相对零量化而始终逐点检查。plan字节记录新格式和QB宽度，因此新解压器仍兼容旧32/64位bitstream。
6. partition之间彼此独立，通过OpenMP并行执行采样、预测、量化、Huffman和Zstd；解压端按相反顺序恢复Huffman符号，并在每个element内依次重建节点。

FDWZ与MRCZ一样在编译期生成正式、分段计时和消融三个版本。正式版本不包含各阶段的计时调用和计时数组；命令行使用`-measure-phase-time true`时，计时版本额外报告sampling以及预测/量化、Huffman和Zstd等阶段的累计partition CPU时间；`FDWZ-ablation`仍使用不计时版本并读取消融开关。

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
| 三维数据中的六种维度顺序 | **采纳为partition级选择** | 不存在跨数据和partition普遍最优的固定顺序，因此由确定性的等间隔sample选择。官方SZ3调优只比较正序与逆序，但FDWZ限制为`XYZ/ZYX`后，即使使用完整数据oracle选择，小数据65个field/EB case的平均BR regret约2.97%、最大约21.80%，暂不照搬。 |
| cubic插值 | **采纳** | 正式预测固定使用cubic；采样阶段比较六种linear顺序并将最佳方向映射到对应cubic，避免为选择方向重复运行昂贵的cubic sample。（SZ3中已有） |
| surface-first：按边、面、体及中心向外顺序预测 | **放弃** | 在代表性数据上没有形成rate--distortion优势。 |
| normal-boundary：先建立中心块，再沿法向外插边界 | **放弃** | 五个数据上均被正式物理中点预测树Pareto支配。 |
| 官方SZ3式2幂stride、单原点seed的index-order | **放弃** | 对小element反复支付粗层级远距离预测成本，三个小数据均无优势。 |
| 不预测物理中点，而是预测索引中点 | **放弃** | pb0和jicf0的BR分别增加约8.4%和4.7%，说明物理距离比索引邻近更重要。 |
| element编号方向预测 | **放弃** | 在五个数据上严格优于所有element内方向的点不足约1.1%，加入空间oracle后的总体改善很小。 |
| 真实共享面邻居与坐标匹配预测 | **放弃当前实现方向** | 早期跨公共面实验没有优于element内预测，同时会引入邻居搜索和拓扑metadata成本。 |
| QoZ式GLL预测路径 | **放弃并删除旧实现** | test0/UZ上的BR明显差于FDWZ，不作为正式优化继续；QoZ对照使用官方API。 |

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
| 32位QB与unpredictable原值逃逸 | **采纳** | `INT32_MIN`作为标记，超出32位范围或浮点重构不安全的点在Huffman流后顺序保存无损原值，不再整块回退64位；正常残差通过relative-EB上界走无逐点判断的快速路径。pb0/P、64核同节点测试中，32位QB在`1e-4`和`1e-5`下比64位实现分别快约11.0%和14.8%；18个小数据field的90个常规case全部通过误差界，极端double逃逸测试也通过。 |
| 从sample直接选择tree数 | **放弃** | 小sample无法正确摊销tree metadata，tree数与完整partition的最优选择经常不一致。 |
| 按符号数阈值或half/full sample外推tree成本 | **放弃** | regret仍不够稳定；最终固定最大有效tree比自适应规则更简单稳健。 |
| 每个partition的多个Huffman流合并到一个Zstd frame | **采纳** | 保留不同符号分布的Huffman收益，同时让Zstd继续利用整个partition的字节级冗余。 |

## 2.6. 配置采样、并行与实现加速阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| partition内确定性的均匀采样 | **采纳interval=128、min=32** | sample比较六种linear顺序，正式压缩映射到同方向cubic。相对完整比较六种固定cubic的oracle，18个小数据field在各自共同PSNR近中点处平均CR regret约0.34%、最大4.34%，均低于5%。 |
| 只采样partition连续前缀 | **放弃** | 对非平稳partition代表性不足，曾出现超过100%的holdout regret。 |
| 预计算prediction stencil和GLL权重 | **采纳** | 避免在每个element重复构造预测关系，并保证编解码路径一致。 |
| partition级OpenMP并行 | **采纳** | partition天然独立，无需在element预测树内部增加同步。 |
| 四个element分批执行正式预测 | **采纳** | 保持每个element内部预测顺序不变，同时复用stencil读取并增加独立指令；两个本地数据上的端到端压缩时间减少约2%--5%，bitstream逐字节不变。八个element一批没有进一步收益。 |
| 未初始化的Huffman/Zstd输出工作区 | **采纳** | 避免`vector(size)`对最坏情况缓冲做无意义清零，仍将实际allocation计入codec wall time；两个本地数据上的端到端压缩时间减少约9%--14%，bitstream逐字节不变。 |
| 跨element AVX2并行预测与量化 | **放弃正式合并** | 约256 KB tile解决了整partition position-major访问的缓存退化，但float没有收益，double端到端仅约快9%，不足以承担独立舍入语义和复杂实现。 |
| 为固定`8x8x8` shape展开直接索引kernel | **放弃** | pb0上保持bitstream和误差结果不变，但总压缩时间只减少约2%--6%，同时生成约3.9万行代码并增大binary，不足以抵消维护成本。 |
| 在同一预测轮内改为逐线遍历 | **放弃** | 相比当前按目标地址递增的memory顺序，总时间变化为约`+0.5%`、`+0.9%`和`-0.3%`，没有稳定收益；保留当前memory顺序。 |
| 用紧凑模板表替代展开的stencil/plan表示 | **放弃** | 总压缩时间慢约0.6%--1.4%，decision阶段慢约1.8%--7.1%，binary也更大。 |

# 3. MRCZ核心与优化记录

MRCZ是一种面向高阶有限元/谱元数据的model-representation-based compressor。它不直接预测某个节点，而是将每个element转换为一组描述局部场的基系数，以量化后的低误差模型重建值作为全部节点的预测，再单独编码逐点residual。这个两级表示同时利用了element内的高阶光滑性，并通过residual量化提供与具体变换无关的严格逐点误差界限。

## 3.1. 核心流程

1. 根据element shape、各维GLL节点和选定基函数，为每个有效维度构造一维基函数求值矩阵及其逆矩阵；正式MRCZ使用定义在真实GLL节点上的DChT（Discrete Chebyshev Transform）。三维正变换和反变换均按三个维度依次执行张量积矩阵运算。
2. 对每个element执行正变换，将nodal values转换为transform coefficients。所有系数使用`coeff_eb=0.05*abs_eb`统一量化，并立即用量化后的系数值覆盖原系数。
3. 对量化系数执行反变换，得到该element所有节点的近似值。随后计算原始值与该近似值之间的pointwise residual，并使用用户指定的`abs_eb`量化。解压端以相同量化系数完成反变换并叠加重建residual，因此最终逐点误差不超过`abs_eb`；变换本身不需要正交，只需要可逆且数值稳定。
4. 每个partition内的coefficient QB由element-major转为coefficient-major，使相同基模式在不同element间连续排列；residual QB保持element-major。两类QB分别建立一棵Huffman tree，Huffman metadata和payload拼接后共同压缩为一个Zstd frame。
5. partition之间独立并通过OpenMP并行。每个线程复用一组partition-local coefficient/residual缓冲区；解压端按Zstd、Huffman、系数反量化、反变换和residual叠加的逆向流程恢复数据。

MRCZ共享同一份模板实现并在编译期生成三个版本：`MRCZ`且`if_measure_phase_time=false`时使用不含分段计时代码的正式版本；`MRCZ`且该开关为`true`时使用带分段计时的版本；`MRCZ-ablation`使用消融版本并读取`param`中的MRCZ优化开关。命令行程序可通过`-measure-phase-time true`选择计时版本。

命令行消融使用`-method MRCZ-ablation`并通过`-grid-assumption`、`-transform`、`-coeff-eb-scale`、`-coefficient-transpose`和`-tiled-coefficient-transpose`设置MRCZ开关。FDWZ消融使用`-method FDWZ-ablation -max-huffman-tree false`切换到单棵Huffman tree，也可用`-fdwz-plan <order>-<linear|cubic>`跳过采样并固定预测方案；正式`MRCZ`和`FDWZ`使用论文采用的固定配置，不受消融开关影响。`-fdwz-sample-interval`不是消融开关，可控制正式FDWZ的采样密度。消融开关不写入metadata，因此解压消融结果时必须提供与压缩相同的参数。FDWZ选中的完整预测方案属于解码必需的算法状态，仍会随该partition的bitstream保存。

## 3.2. 变换与coefficient error bound阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| 在真实GLL节点上使用DChT | **采纳** | 旧稿实验中，同PSNR下DChT+GLL相对DLT+GLL约改善20% CR；主要收益来自更集中的系数分布，而不是residual差异。 |
| DLT、DST、DCT以及将节点视为uniform的变换 | **放弃为正式配置** | 已比较不同基函数和节点假设，整体rate--distortion不如DChT+GLL稳定；实现仍允许构造部分变换，但正式MRCZ固定使用DChT+GLL。 |
| 全局`coeff_eb_scale=0.05` | **采纳** | 在系数码率与反变换预测精度之间形成稳定折中；过小会浪费用户误差预算，过大会增加BR。 |
| 按模式的L1/RMS敏感度、乘积或总阶数分配coefficient EB | **放弃** | 不同数据上的BR变化方向不一致，系数流与residual流经常相互抵消，未形成稳定规则。 |
| 扩大量化为零的低频系数区间 | **放弃** | mixlay0和jicf0上的同PSNR结果均更差，不能仅依靠增加零系数改善完整bitstream。 |

## 3.3. coefficient与pointwise residual阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| partition级coefficient zero map，用于标注某个model的系数是否在partition内所有element上全部被量化为0 | **放弃** | jicf0五个field、五个EB上的CR平均仅提高约0.107%、最大约0.486%，不足以补偿扫描、metadata和格式复杂度。 |

## 3.4. QB布局、Huffman与Zstd阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| coefficient-major布局 | **采纳** | 相同基模式跨element连续排列，可聚集高阶模式中的大量零值及同模式统计规律；旧稿实验中约改善20% CR。 |
| coefficient与residual各使用一棵Huffman tree | **采纳** | 两类QB分布不同，独立建树后仍可放入同一Zstd frame利用partition级字节冗余。 |
| 按系数总阶数拆成2--4棵Huffman tree | **放弃** | 小数据上最多有约10%的平均BR收益，但jicf0大partition上三树和四树相对单树仅改善约1.1%且依赖数据区域，正式格式保持单棵coefficient tree。 |
| coefficient和residual的Huffman流合并到一个Zstd frame | **采纳** | 避免为同一partition支付多个Zstd frame，并允许Zstd继续利用两段Huffman数据的字节级冗余。 |

## 3.5. 并行、内存与实现加速阶段

| 优化或候选方案 | 状态 | 结论 |
|---|---|---|
| 支持任意`(n0,n1,n2)`的统一变换内核 | **采纳** | 等尺寸专用内核在8x8x8 element上的性能差异不足1%，删除重复的regular/irregular路径更易维护。 |
| 连续指针访问并由编译器自动生成AVX2 | **采纳** | 重构六个张量变换内层循环后，GCC能够使用32-byte vector；pb0/T计算节点样例的压缩和解压分别加速约1.36x和1.38x，bitstream与误差结果不变。 |
| 32x32分块coefficient transpose | **采纳** | 以cache-aware tiled transpose完成element-major到coefficient-major的布局转换，不改变bitstream；pb0/T样例中压缩和解压时间分别减少约11.3%和14.7%。IPDPS路线可将其作为数据布局与缓存优化的算法贡献，DCC路线只作为实现加速，不单列核心算法贡献。 |
| 只添加`#pragma omp simd`或手写transform intrinsics | **放弃** | 单独添加pragma曾使压缩变慢约1.5%--3%；指针重构已能触发稳定自动向量化，无需维护手写intrinsics。 |
| 跨element batched transform | **放弃正式合并** | float/double批处理能够逐字节复现正式bitstream，但独立消融中transform本身反而变慢；组合收益主要来自绕过素朴转置，已由32x32分块转置以更简单方式获得。 |
| partition级OpenMP与独立`partition_num` | **采纳** | partition边界决定bitstream，`max_cpus`只决定并发度；固定partition时不同线程数可产生相同结果，并避免实现partition内并行Huffman。 |
| 每线程复用partition-local QB缓冲区 | **采纳** | 避免同时为整个field保留两份完整`int64_t`数组；pb0/T样例的峰值RSS由约595 MB降至245 MB，同时略微缩短压缩时间。 |

# 4. 消融实验

- 调优coefficient transpose的tile size：比较`8x8`、`16x16`、`32x32`、`64x64`和`128x128`，分别记录压缩与解压transpose阶段时间；当前`32x32`只验证过优于未分块transpose，尚未证明是最优值。同步检查各版本是否由编译器生成SIMD/AVX2指令，区分cache blocking与自动向量化的收益。
- 测试SZ2的`SZ_compress_customize_threadsafe`：在进入OpenMP并行区前完成一次`SZ_Init()`，为各partition使用相同数据类型和ABS EB并恢复partition级并行，解压端暂时保持串行。使用不同局部值域的partition反复检查bitstream可解码性、逐点误差界限和结果稳定性，通过后再用于正式throughput实验。

# 5. Rate-distortion结果

以下CR增长均由RD曲线在相同PSNR下插值得到，不比较相同REL EB的点。对每个field先求MRCZ、FDWZ和五个正式官方baseline共同覆盖的PSNR区间，再选择最接近区间中点的整数PSNR；不同field使用各自的目标PSNR。CR增长由BR之比计算，例如MRCZ相对FDWZ为`BR_FDWZ/BR_MRCZ-1`。pb0、jicf0和小数据的FDWZ数字均已于2026-09-05使用当前linear-proxy、`fdwz_sample_interval=128`和32位QB/unpredictable版本重跑。

| 数据 | 比较 | CR平均增长 | CR最大增长 | 最大值对应field |
|---|---|---:|---:|---|
| pb0（当前FDWZ） | MRCZ相对FDWZ | 125.30% | 130.00% | P |
| pb0（当前FDWZ） | FDWZ相对最优官方baseline | 37.84% | 49.05% | UX |
| jicf0（当前FDWZ） | MRCZ相对FDWZ | 97.49% | 131.83% | P |
| jicf0（当前FDWZ） | FDWZ相对最优官方baseline | 123.72% | 205.87% | P |
| 大数据（当前FDWZ） | MRCZ相对FDWZ | 111.39% | 131.83% | jicf0/P |
| 大数据（当前FDWZ） | FDWZ相对最优官方baseline | 80.78% | 205.87% | jicf0/P |
| 小数据（当前FDWZ） | MRCZ相对FDWZ | 136.66% | 669.47% | tc_channel0/UY |
| 小数据（当前FDWZ） | FDWZ相对最优官方baseline | 58.83% | 359.95% | MFEM ex0/no_name_field |

五个baseline为SZ3-4D、SZ2-3D-fused、QoZ-4D、ZFP-3D-fused和SPERR-2D；每个field在目标PSNR处取其中最低BR。大数据汇总10个field；小数据不计入RD曲线接近PSNR轴的test0/UZ，汇总其余17个field。各field等权，逐field共同区间、目标PSNR和最优baseline保存在`results/RD/midpoint_cr_summary.csv`，统计可由`python3 py_funcs/summarize_rd.py`复现。

## 5.1. 小数据布局选择

小数据实验使用mixlay0、test0、tc_channel0三个NekRS数据集、MFEM ex0/ex1两个样例和Nektar++ cylinder样例。各数据的非常数field如下：

| 数据 | 非常数field |
|---|---|
| mixlay0 | UX、UY、P、T、S00、S01 |
| test0 | UZ、T、S00、S01 |
| tc_channel0 | UX、UY、UZ |
| MFEM ex0 | no_name_field |
| MFEM ex1 | no_name_field |
| Nektar++ cylinder | p、u、v |

实验比较MRCZ、FDWZ以及官方SZ3、SZ2、QoZ、ZFP和SPERR API。每个case使用完整field、1 CPU和1个partition，REL EB为`1e-2`、`1e-3`、`1e-4`、`1e-5`和`1e-6`；正式结果重测同样固定使用1 CPU。各官方API可用的布局如下：

| 官方API | 1D | 2D | fused 3D | per-element 3D | 4D |
|---|:---:|:---:|:---:|:---:|:---:|
| SZ3 | ✓ | ✓ | ✓ | ✓ | ✓ |
| SZ2 | ✓ | ✓ | ✓ | ✓ | ✓ |
| QoZ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ZFP | ✓ | ✓ | ✓ | ✓ | ✓ |
| SPERR | — | ✓ | ✓ | ✓ | — |

设一个partition包含`E`个element，每个element的shape为`[N0,N1,N2]`，`N=N0*N1*N2`。各布局调用官方API的方式如下：

| 布局 | API调用 |
|---|---|
| 1D | 将整个partition作为`[E*N]`调用一次 |
| 2D | 作为`[E,N]`调用一次 |
| fused 3D | 作为`[E*N0,N1,N2]`调用一次 |
| per-element 3D | 对每个`[N0,N1,N2]` element分别调用，并保存各bitstream长度 |
| 4D | 作为`[E,N0,N1,N2]`调用一次 |

官方库中需要保留的特殊处理如下；SPERR缺少1D和4D入口已由上表体现。

| 官方API | 特殊处理 |
|---|---|
| SZ3 | float路径传入`0.999*abs_eb`，避免曾观测到的约0.031%误差越界 |
| SZ2 | 会删除长度为1的维度；完整3D element的4D与fused 3D路径产生相同bitstream。其状态为进程级全局变量，不能安全并发，因此正式封装串行调用且不报告多CPU scaling |
| SZ2（MFEM ex1） | `5x5x1` per-element布局在REL `1e-4`会触发官方库allocator断言；失败保留在原始结果中 |
| QoZ | 删除长度为1的维度；q4自动调优只用于2D/3D，1D通过官方API固定使用插值后端，以避开默认路径的heap corruption |

完整布局矩阵包含25种方法、18个非常数field和5个REL EB，共2250个case；其中2249个成功并满足error bound，唯一失败是上表所列MFEM ex1的per-element SZ2案例。

在各布局共同覆盖的PSNR范围内插值BR并比较regret后，正式大数据baseline固定为SZ3-4D、SZ2-3D-fused、QoZ-4D、ZFP-3D-fused和SPERR-2D。per-element 3D通常因大量独立bitstream、长度metadata及无法利用跨element冗余而明显退化；test0/UZ的RD曲线接近PSNR轴且无法有效区分布局，因此保留原始结果但不用于选择。这些布局在pb0和jicf0实验前即已冻结，不使用大数据结果反向挑选baseline。

## 5.2. 小数据正式方法结果

下图展示全部18个非常数field上25种方法的BR--PSNR曲线；同一种官方压缩器的不同API布局使用相同颜色和不同线型，BR横轴使用线性尺度并固定从0开始。各子图聚焦主要曲线的有效比较范围，允许高BR方法和明显偏高的PSNR端点在坐标边界处截断。FDWZ曲线于2026-09-05使用当前cubic-only、linear-proxy、`fdwz_sample_interval=128`、32位QB/unpredictable和1 CPU/1 partition重跑；90个case全部满足误差界且没有触发unpredictable，其余24种方法沿用已冻结的正式结果。

![Small-data rate-distortion](results/RD/small_rd.svg)

## 5.3. 大数据正式结果

大数据实验使用pb0的UX、UY、UZ、P和T五个field，以及jicf0的UX、UY、UZ、P和S00五个field。正式方法为MRCZ、FDWZ、SZ3-4D、SZ2-3D-fused、QoZ-4D、ZFP-3D-fused和SPERR-2D；pb0使用64 CPU和64个partition，REL EB为`1e-2`、`1e-3`、`1e-4`、`1e-5`和`1e-6`；jicf0使用128 CPU和128个partition，REL EB为`1e-4`、`1e-5`、`1e-6`、`1e-7`和`1e-8`。

pb0的FDWZ于2026-09-05在64 CPU、64个partition的独占节点上完成当前版本重测（Slurm job `20422881`）；五个field、五个REL EB共25个case均通过项目的`100.01%`浮点数值容差检查（最大为`100.003336%`）且没有触发unpredictable。jicf0于同日在128 CPU、128个partition的独占节点上重测；主任务`20423381`在完成前四个REL EB后达到45分钟时限，续跑任务`20424117`补齐`1e-8`。jicf0的25个case同样全部通过容差检查（最大`100.000000%`）且没有触发unpredictable。按照本章开头的口径，MRCZ相对FDWZ的CR平均增长在pb0和jicf0上分别为125.30%和97.49%；FDWZ相对最优官方baseline的CR平均增长分别为37.84%和123.72%。

pb0的当前64 CPU结果如下；BR横轴使用线性尺度并从0开始，允许明显偏高的曲线端点在坐标边界处截断。MRCZ和五个官方baseline沿用冻结结果，仅FDWZ使用当前版本重跑。

![pb0 rate-distortion](results/RD/pb0_rd.svg)

jicf0的当前128 CPU结果如下；其他六种方法同样沿用已冻结结果。

![jicf0 rate-distortion](results/RD/jicf0_rd.svg)

# 6. Throughput

pb0和jicf0的正式throughput分别使用64 CPU和128 CPU测量。吞吐量均为原始数据量除以backend内部耗时，并按CPU数归一化为MB/s per CPU，纵轴使用对数尺度；backend为压缩或解压API调用本身，其内部QB、Huffman和Zstd工作区allocation仍计入对应时间。两个数据的FDWZ使用上述当前版本结果，其他六种方法沿用冻结结果。pb0的25个同REL EB点上，FDWZ压缩和解压吞吐均全部高于SZ3-4D，几何平均吞吐比分别为1.333和1.396，最小比分别为1.083和1.278。jicf0上，FDWZ压缩在25点中的18点快于SZ3-4D，几何平均吞吐比为1.059；解压在25点中的24点更快，几何平均吞吐比为1.482。因此FDWZ的解压优势在两个大数据上都很稳定，压缩端则只能说总体达到SZ3-4D速度，jicf0的严格误差界限下仍会明显落后。另在独占节点上对pb0/P的`1e-4`与`1e-5`各重复三次，FDWZ压缩分别比SZ3-4D快19.8%和20.1%，解压分别快30.7%和37.6%。

pb0结果首先以实验直接指定的REL EB为横坐标；同一横坐标位置上的所有方法具有相同的error bound和计算精度。

![pb0 compression throughput by REL EB](results/RD/pb0_compression_throughput_eb.svg)

![pb0 decompression throughput by REL EB](results/RD/pb0_decompression_throughput_eb.svg)

随后以各方法实际达到的PSNR为横坐标，从fixed-quality角度展示压缩与解压throughput。

![pb0 compression throughput](results/RD/pb0_compression_throughput.svg)

![pb0 decompression throughput](results/RD/pb0_decompression_throughput.svg)

jicf0以REL EB和实际PSNR为横坐标的对应结果如下。

![jicf0 compression throughput by REL EB](results/RD/jicf0_compression_throughput_eb.svg)

![jicf0 decompression throughput by REL EB](results/RD/jicf0_decompression_throughput_eb.svg)

![jicf0 compression throughput](results/RD/jicf0_compression_throughput.svg)

![jicf0 decompression throughput](results/RD/jicf0_decompression_throughput.svg)

# 7. Weak scaling

实验将pb0和jicf0分别均分为64和128个全局partition，每个CPU固定处理一个partition。pb0测试`1、2、4、8、16、32、64` CPU，jicf0测试`1、2、4、8、16、32、64、128` CPU。N个CPU处理连续的N个partition，因此CPU数和数据量同步增加。pb0每个partition约含30424个element和59.42 MiB原始数据，jicf0约含25315个element和96.56 MiB原始数据；目标PSNR分别为100 dB和150 dB。

CPU数量较少时，实验只处理完整field中的局部区间，而不同区间的data range可能不同。如果直接使用REL EB，程序会根据局部data range得到不同的ABS EB，使各CPU点的实际误差界和重建质量不再一致。因此，每个method/field先根据完整数据的RD曲线插值得到目标PSNR对应的REL EB，再乘以完整field的data range得到ABS EB；所有CPU数量和抽样区间均使用这个固定的ABS EB。同时，实验通过`-compute-precision`固定使用完整数据在该REL EB下选择的float或double，避免局部data range变化导致计算精度切换。局部区间的data range还可能小到触发constant compression，因此weak scaling压缩命令统一设置`-constant-compression false`，确保每个case都实际运行对应backend。汇总PSNR也使用完整field的data range和各区间的加权MSE重新计算。

每个CPU数、每个field测试4次。先将完整数据固定划分为最大CPU数个全局partition；使用N个CPU时，一个测试段由N个相邻partition组成，再从所有可选测试段中等间隔选择4段。例如jicf0使用2 CPU时选择全局partition `[0,2)`、`[42,44)`、`[84,86)`和`[126,128)`。当N较大、完整数据只能提供2段或1段时，每段分别重复2次或4次。因此每个CPU点仍覆盖jicf0中多个位置，而不是复制某个代表性partition。每次测试先汇总全部五个field的原始数据量和内部时间，再计算综合throughput与efficiency。图中展示四次测试合并后的综合efficiency，不绘制误差棒；四次测试的样本标准差仍保存在summary CSV中。Slurm脚本使用`OMP_PLACES=cores`和`OMP_PROC_BIND=spread`将OpenMP线程固定到物理核，并使用`--exclusive --mem=235G`独占节点及235 GiB内存，避免并发Huffman工作区触发cgroup OOM或受到同节点任务干扰；pb0程序仍只使用64个线程，jicf0使用128个线程。

图中各方法的efficiency以自身1 CPU吞吐量为基准；SZ2官方实现使用进程级全局状态，不能安全并发，因此只有1 CPU结果且不画入趋势图。pb0的FDWZ于2026-09-05在节点a173上使用当前版本重跑（Slurm job `20424829`），140个raw case均通过`100.01%`误差容差，各field/CPU汇总PSNR为99.942--100.015 dB。其他方法沿用节点a214上统一线程绑定的冻结结果；efficiency均由同一方法在同一节点上的结果自归一化，因此节点间绝对速度差异不会直接进入efficiency。FDWZ的压缩/解压efficiency在8 CPU时为0.880/0.841，64 CPU时为0.241/0.227；由于其单核速度更高，它比其他方法更早撞上内存带宽瓶颈。即便如此，64 CPU绝对压缩与解压吞吐仍为4422和5523 MB/s，分别比冻结的SZ3-4D结果高18.7%和38.1%；该幅度大于已观察到的约2%--5%跨节点差异。jicf0的当前FDWZ结果尚待重测。

![pb0 weak-scaling efficiency](results/weak_scaling/weak_pb0_efficiency.svg)

# 8. Strong scaling

# 9. End-to-end test

# 10. 当前状态

- 正式MRCZ算法已冻结：逐element DChT、系数与残差量化、32x32分块coefficient-major
  转置、Huffman和Zstd，partition之间OpenMP并行。
- FDWZ使用最大有效Huffman tree数和cubic正式预测，并在每个partition以128个element为默认采样间隔、至少32个均匀样本比较六种linear顺序，再映射到同方向cubic；QB固定为32位，超界或重构不安全的值通过unpredictable原值流逃逸。
- 官方API baseline及tensorization已经冻结，定义见第5.1节。
- README中的18个小数据field FDWZ曲线已使用当前linear-proxy与32位QB/unpredictable实现重跑，90个case全部满足误差界；其余24种方法沿用冻结结果。
- pb0与jicf0的rate--distortion与throughput均已使用当前FDWZ版本重跑，各25个case都通过`100.01%`浮点数值容差检查且没有触发unpredictable；其他六种方法沿用冻结结果。
- pb0 weak scaling已使用当前FDWZ版本重跑，140个case均通过误差容差；jicf0的当前FDWZ结果待重测。

# 11. 投稿前主线

[ ] 1. 冻结实验版本与结果。
       提交当前正式代码、实验脚本和结果说明；记录生成正式CSV所用的commit、Slurm job、
       编译器、CPU、partition数量和计时口径。CSV继续保存在对应数据旁边。

[ ] 2. 生成论文用rate-distortion汇总。
       从pb0和jicf0各field CSV生成BR--PSNR图，并按第5章的逐field共同区间近中点口径
       汇总MRCZ、FDWZ和五个官方baseline的BR、CR、最大误差与异常点。旧稿中的数字不得直接复用。

[ ] 3. 整理速度与扩展性结果。
       绘制single-core throughput、fixed-quality strong/weak scaling、speedup和parallel
       efficiency；明确64--128 CPU后的内存带宽与NUMA瓶颈，不声称近线性扩展。
       所有正式性能实验全局统一使用物理核与OpenMP线程绑定：Slurm申请
       `--hint=nomultithread`，`srun`使用`--cpu-bind=cores`，并设置
       `OMP_DYNAMIC=false`、`OMP_PLACES=cores`和`OMP_PROC_BIND=spread`；记录这些设置，
       禁止混用未绑定线程的旧throughput或scaling结果。
       Strong scaling的所有CPU点必须固定element范围、error bound和partition网格，在同一
       Slurm allocation中运行，并核验申请的CPU数与程序实际使用的CPU数一致。

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

# 12. 投稿后或artifact阶段

[ ] 7. 审计并版本化bitstream metadata，增加未知版本、截断文件、错误method和错误数据类型测试。消融配置不写入metadata，由用户保证编解码参数一致。

[ ] 8. 让t1在超过error bound、重建数量错误或编解码失败时自动返回非零状态，并处理低于机器精度的EB。

[ ] 9. 清理正式路由中剩余的实验方法和历史注释。

[ ] 10. 为未显式指定`partition_num`时定义与`max_cpus`无关的默认规则，并在jicf0上验证固定partition时不同线程数产生相同bitstream、BR、PSNR和最大误差。

# 13. 不进入当前论文

- 坐标X/Y/Z专用有损压缩及几何/拓扑安全定义。
- NekRS尾metadata的语义解析；当前继续作为opaque blob无损保存。
- coefficient delta、更多SZ3预测方案、跨element AVX2和GPU实现。
