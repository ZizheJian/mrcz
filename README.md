# MRCZ

MRCZ研究高阶有限元与谱元模拟数据的误差有界有损压缩，主要面向NekRS/Nek5000等程序生成的三维物理场，并使用MFEM和Nektar++样例验证数据布局的通用性。

此类数据由大量高阶element组成。element之间按非规则网格及拓扑关系组织，每个element内部的采样点位于Gauss--Lobatto--Legendre（GLL）节点上，形成尺寸较小但结构固定的张量积网格。项目包含两个针对这种数据组织设计的压缩器：

- FDWZ（Fast Distance-Weighted Interpolation-Based Lossy Compressor）在每个element内按照GLL节点距离执行分层插值预测。
- MRCZ（Model-Representation-Based Compressor）将每个element变换为局部模型系数，并通过逐点residual correction提供误差界限。

二者均以partition为独立压缩单元，通过OpenMP并行处理partition，并在一个partition内共享Huffman和Zstd编码。

## 安装

需要Python 3、CMake、支持C++20和OpenMP的C++编译器。绘图与分析脚本还需要NumPy、SciPy、Matplotlib、Pillow和Plotly。

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
python -m pip install numpy scipy matplotlib pillow plotly
```

MGARD-X还需要Protobuf。在Anvil上构建前加载对应模块：

```bash
module load protobuf/3.11.4
```

主实验入口为`t1_start_c_code.py`。它读取JSON配置，构建所需依赖并运行压缩与解压实验：

```bash
python3 t1_start_c_code.py CONFIG.json
```

如果只需要构建，可以手动运行：

```bash
cmake -S third_party/sperr -B build/sperr_build -DCMAKE_INSTALL_PREFIX="$PWD/install/sperr_install" -DBUILD_SHARED_LIBS=OFF
cmake --build build/sperr_build --parallel 8
cmake --install build/sperr_build

cmake -S third_party/mfem -B build/mfem_build -DCMAKE_INSTALL_PREFIX="$PWD/install/mfem_install"
cmake --build build/mfem_build --parallel 8
cmake --install build/mfem_build

cmake -S third_party/MGARD -B build/mgard_serial -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install/mgard_install" -DMGARD_ENABLE_SERIAL=ON -DMGARD_ENABLE_OPENMP=OFF -DMGARD_ENABLE_CUDA=OFF -DMGARD_ENABLE_HIP=OFF -DMGARD_ENABLE_SYCL=OFF -DMGARD_ENABLE_MDR=OFF -DMGARD_ENABLE_UNSTRUCTURED=OFF -DMGARD_ENABLE_BENCHMARKS=OFF -DMGARD_ENABLE_DOCS=OFF -DBUILD_TESTING=OFF
cmake --build build/mgard_serial --parallel 8
cmake --install build/mgard_serial

cmake -S c_codes -B build -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel 8
cmake --install build
```

安装后的命令行程序位于`install/bin/mrcz_compress`和`install/bin/mrcz_decompress`。外部C++项目也可以通过CMake导入header-only核心：

```cmake
find_package(MRCZ CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE mrcz::mrcz)
```

公共header为：

```cpp
#include <mrcz_core/mrcz.hpp>
```

## FDWZ

FDWZ借鉴分层插值预测，但将数值依赖限制在单个element内，并使用reference-space GLL坐标构造预测顺序和插值权重。

当前流程如下：

1. 根据element shape与各维GLL坐标预先构造预测层、访问顺序、reference位置和插值权重。每个有效维度的两个端点构成初始边界，随后递归选择reference space中最接近区间中点的未编码节点。
2. 二维或三维element的4或8个角点以零为预测值直接量化。其余节点使用当前可用的重建reference执行三次、二次或一次Lagrange插值。
3. 每次预测后立即量化residual并覆盖目标值，后续预测只引用重建值，使压缩端与解压端采用相同的误差传播路径。
4. 每个partition以固定间隔均匀采样element，比较六种linear维度顺序，再以选中的方向执行正式cubic预测。默认采样间隔为128，最少采样32个element。
5. 量化符号按照角点和预测层分组，并使用element shape允许的最大有效Huffman tree数量，最多五棵。所有Huffman metadata与payload合并到一个Zstd frame。
6. 量化符号使用32位整数。`INT32_MIN`标记无法安全表示的值，对应原始浮点值按出现顺序附加到Huffman流之后。
7. partition之间独立并通过OpenMP并行。预测stencil与GLL权重预先计算，正式预测每次批量处理四个element。

## MRCZ

MRCZ以每个element的局部变换表示作为预测模型，再单独编码逐点residual。当前流程如下：

1. 根据element shape和各维GLL节点构造DChT（Discrete Chebyshev Transform）的一维基函数求值矩阵及逆矩阵。
2. 对每个element执行张量积正变换。所有系数使用`coeff_eb=0.05*abs_eb`量化，并立即以量化后的值覆盖原系数。
3. 对量化系数执行逆变换，将重建模型作为全部节点的预测值。原始值与预测值之间的pointwise residual使用用户指定的`abs_eb`量化。
4. coefficient QB从element-major重排为coefficient-major，使相同变换模式在不同element间连续排列。residual QB保持element-major。
5. coefficient与residual分别建立一棵Huffman tree，两个流共同压缩为一个Zstd frame。
6. 全局coefficient-major重排通过最多覆盖32个element和32个mode的cache-aware tile完成。每个OpenMP线程复用partition-local coefficient与residual缓冲区。

## 数据与baseline布局

仓库submodule包含三个可直接测试的NekRS小数据：

- `third_party/pySEMTools/examples/data/mixlay0.f00001`
- `third_party/pySEMTools/examples/data/tc_channel0.f00001`
- `third_party/fld_pyutils/demos/data/test0.f00001`

MFEM和Nektar++验证数据可重新生成：

```bash
python3 generate_mfem_and_nektar_data.py
```

脚本生成两个MFEM样例及Nektar++ cylinder VTU数据。论文中的两个大规模NekRS数据集不保存在Git仓库中，其位置通过实验JSON配置指定。

baseline通过未经修改的官方API接入。设一个partition包含`E`个element，每个element shape为`[N0,N1,N2]`，且`N=N0*N1*N2`，候选布局定义如下：

| 布局 | API调用方式 |
|---|---|
| 1D | 将整个partition作为`[E*N]`调用一次 |
| 2D | 作为`[E,N]`调用一次 |
| fused 3D | 作为`[E*N0,N1,N2]`调用一次 |
| per-element 3D | 对每个`[N0,N1,N2]` element分别调用，并保存各bitstream长度 |
| 4D | 作为`[E,N0,N1,N2]`调用一次 |
| 4D+GLL | MGARD-X的4D nonuniform-coordinate接口 |

论文中的正式大数据baseline布局为SZ3-4D、SZ2-3D-fused、QoZ-4D、ZFP-3D-fused、SPERR-2D和MGARD-X-1D。SZ2的官方API使用进程级可变状态，不能安全执行并发调用，因此不参与多核throughput实验。MGARD-X的内存需求不允许其在大规模数据上运行与其他方法相同的高核数配置，相关图中只保留实测的小核结果和理想线性扩展上界。

## 实验与绘图

JSON配置保存在`configs/`，实验CSV和绘图中间结果保存在`results/`。主要入口包括：

| 文件 | 用途 |
|---|---|
| `t1_start_c_code.py` | 构建并运行rate-distortion、throughput和重建实验 |
| `t2_find_three_irregular_elements.py` | 在NekRS网格中寻找三个两两共享面的非平面element，并生成用于论文几何示意图的Plotly HTML |
| `t3_test_cubic_face_prediction.py` | 比较共享面上的跨element与面内cubic预测，汇总误差和量化符号熵；由prediction-entropy Slurm脚本调用 |
| `generate_mfem_and_nektar_data.py` | 生成MFEM与Nektar++验证数据 |
| `py_funcs/summarize_rd.py` | 在共同PSNR区间汇总逐field rate-distortion结果 |
| `py_funcs/summarize_end_to_end.py` | 合并codec时间与Globus传输记录 |
| `py_funcs/plot_style.py` | 统一论文图的字体、方法顺序、配色、marker和坐标轴样式 |
| `py_funcs/rasterize_paper_figures.py` | 按最终物理尺寸将复杂源图转换为300 dpi PNG |

主要论文图可由以下命令重新生成：

```bash
python3 py_funcs/plot_fdwz_tree_ablation.py
python3 py_funcs/plot_fig8_transform_comparison.py
python3 py_funcs/plot_transform_mode_statistics.py
python3 py_funcs/plot_mrcz_reordering_ablation.py
python3 py_funcs/plot_small_rd.py
python3 py_funcs/plot_large_rd.py
python3 py_funcs/plot_combined_throughput.py
python3 py_funcs/plot_end_to_end.py
python3 py_funcs/plot_transfer_size_time.py
```

可视化图还需要重建场和ParaView切片中间文件：

```bash
reconstruct_job=$(sbatch --parsable run_visualization_cr33_reconstruct.sbatch)
sbatch --dependency=afterok:${reconstruct_job} run_visualization_cr33_render.sbatch
```

渲染任务会在`results/visualization/cr33_panel_masters/`保存原图和七种方法的1200 dpi单子图母版，并从这些母版缩小生成300 dpi论文图。以后只调整版式时可以跳过重建和ParaView切片：

```bash
python3 py_funcs/plot_pb0_cr33_visualization.py --compose-only
```

手工流程图直接使用定稿PNG，`cross-element.pdf`可按最终单栏宽度重新生成300 dpi PNG。Matplotlib统计图和场可视化按论文宽度和8 pt字号直接生成300 dpi PNG：

```bash
python3 py_funcs/rasterize_paper_figures.py
```

论文源码位于`overleaf/`，它是Overleaf项目的git克隆，不纳入本仓库。在其中运行`make`可用texlive/20200406编译出`overleaf/main.pdf`。

最终性能数字、消融结果、实验平台和统计口径以论文为准，README不再重复维护这些内容。
