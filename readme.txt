MRCZ
====

克隆项目后，先下载submodule：

    git submodule update --init --recursive

third_party中的pySEMTools和fld_pyutils提供了三个可直接用于MRCZ测试的NekRS数据：

- third_party/pySEMTools/examples/data/mixlay0.f00001
- third_party/pySEMTools/examples/data/tc_channel0.f00001
- third_party/fld_pyutils/demos/data/test0.f00001

t1_start_c_code.py中列出了这三个数据对应的field和error-bound测试配置。当前启用
mixlay0.f00001，其余两个样例可通过取消对应配置的注释来使用。
