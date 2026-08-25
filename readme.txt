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
run_t1.sbatch提供Anvil上的128 CPU批处理入口。

当前开发计划见TODO.txt。MRCZ正式实验保存在experiments/mrcz/，已暂停的SZ3-3D-GLL设计空间
实验及其结论保存在experiments/sz3_3d_gll/。
