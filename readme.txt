MRCZ
====

克隆项目后，先下载submodule：

    git submodule update --init --recursive

third_party中的pySEMTools和fld_pyutils提供了三个可直接用于MRCZ测试的NekRS数据：

- third_party/pySEMTools/examples/data/mixlay0.f00001
- third_party/pySEMTools/examples/data/tc_channel0.f00001
- third_party/fld_pyutils/demos/data/test0.f00001

t1_start_c_code.py中列出了这三个数据对应的field和error-bound测试配置。小数据配置目前作为
注释模板保留；正式大数据通过MRCZ_T1_DATASET选择pb0或jicf0，通过MRCZ_T1_METHODS选择方法。
run_t1_pb0_rd.sbatch和run_t1_jicf_rd.sbatch分别提供Anvil上两个大数据的128 CPU批处理入口。

当前开发计划以及FDWZ设计空间的核心结论见TODO.md。MRCZ正式实验保存在experiments/mrcz/。
