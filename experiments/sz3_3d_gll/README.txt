SZ3-3D-GLL experiments
======================

# XXX 这些文件用于保留2026-08-21到2026-08-24期间的实验过程，不参与正式构建。
# XXX 从项目根目录运行Python脚本；脚本会把新CSV或图片写回本目录。
# XXX 一些历史脚本中的SZ3-3D-GLL和SZ3-3D-GLL-experimental代表当时的2/3/4/5-tree版本，
# XXX 不等于当前已经重新冻结的1tree SZ3-3D-GLL和max-tree SZ3-3D-GLL-optimized。

文件分组
========

- e10：角点可预测性和角点表示方法。
- e11：非角点量化残差按预测访问顺序保存。
- e12--e19：1到5棵Huffman tree的真实bitstream、rate-distortion和大数据抽样。
- e20：六种固定维度顺序。
- e21--e24：GLL距离加权的三次/边界二次插值。
- e25：element内位置、最佳预测方向及element编号方向的可预测性和可视化。
- e26：顺序、linear/cubic和tree数的partition采样选择。
- e27：tree选择regret；最终决定不采样tree数，固定使用shape允许的最大tree数。
- e28：每个局部位置的residual bias；单bin峰值和稳健宽窗口两版均被淘汰。
- run_e*.sbatch和slurm-mrcz-e*：历史计算节点作业及输出。
- OFFICIAL_SZ3_COMPARISON.txt：当前baseline与third_party/SZ3官方实现的语义差异。

最终结论
========

- baseline使用1tree，并按每个element的实际预测访问顺序写入角点和非角点量化值。
- optimized固定使用shape允许的最大tree数，并在每个partition采样选择维度顺序和linear/cubic预测；
  tree数不参与采样。
- vectorization_report.txt记录GCC自动向量化报告；run_avx2_large.sbatch和对应Slurm输出记录跨element
  AVX2原型的float/double测试。约256KB tile消除了整partition position-major访问造成的缓存退化，
  但float没有收益，double端到端约快9%。
- 不继续角点变换、跨element residual转置、特殊prediction tree、reference小EB、
  element编号预测、真实共享面预测和position residual bias。
- 跨element AVX2的最终tile版结果见slurm-sz3-avx2-20135330.out；不合并到正式optimized方法。
- 从2026-08-25起暂停SZ3-3D-GLL的大规模对照、消融、int32 QB存储和AVX2后续优化。
  DCC论文只以MRCZ为主要算法；若以后恢复SZ3研究，从本目录的脚本和结果继续。
