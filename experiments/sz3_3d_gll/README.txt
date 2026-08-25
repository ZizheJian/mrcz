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
- optimized固定使用shape允许的最大tree数，tree数不采样。
- 维度顺序和linear/cubic的采样具有潜在收益，但尚未纳入正式optimized方法。
- 不继续角点变换、跨element residual转置、特殊prediction tree、reference小EB、
  element编号预测、真实共享面预测和position residual bias。
- 后续只需进行正式optimized实现、大规模对照和逐项消融；当前计划见项目根目录TODO.txt。
