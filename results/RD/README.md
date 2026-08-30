# t1 rate--distortion与throughput结果

本目录保存`t1_start_c_code.py`生成的CSV在MRCZ主仓库中的版本化快照。程序仍将运行时结果写在原始数据旁边；完成一轮实验并检查所有方法、error bound和field后，再把完整CSV同步到这里。

- `pb0/`：源自`/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/`中的`pb0.f00001_<field>.csv`和`pb0.f00001_<field>_throughput.csv`。
- `jicf0/`：源自`/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/`中的`jicf0.f00000_<field>.csv`和`jicf0.f00000_<field>_throughput.csv`。
- `test0/`：源自`third_party/fld_pyutils/demos/data/`中的`test0.f00001_<field>.csv`和`test0.f00001_<field>_throughput.csv`。
- `mixlay0/`：源自`third_party/pySEMTools/examples/data/`中的`mixlay0.f00001_<field>.csv`和`mixlay0.f00001_<field>_throughput.csv`。
- `tc_channel0/`：源自`third_party/pySEMTools/examples/data/`中的`tc_channel0.f00001_<field>.csv`和`tc_channel0.f00001_<field>_throughput.csv`。
- `mfem_ex0/`：源自`build/mfem_build/examples/ex0examples/`中的`sol.gf_no_name_field.csv`和对应throughput CSV。
- `mfem_ex1/`：源自`build/mfem_build/examples/ex1examples/`中的`sol.gf_no_name_field.csv`和对应throughput CSV。
- `nektar_cylinder/`：源自`third_party/flow-stability-cylinder/tutorial/base/`中的`Cylinder-Base.vtu_<field>.csv`和对应throughput CSV。

pb0和jicf0是正式大数据结果，各包含五个field，每个field分别保存rate--distortion和throughput，共十个CSV。test0、mixlay0、tc_channel0、mfem_ex0、mfem_ex1和nektar_cylinder是小数据布局选择与实现验证结果。2026-08-30首次复制到本目录；这次复制只建立结果快照，不表示生成这些CSV时使用的commit已经完成审计。

小数据布局实验的汇总结论保存在`experiments/official_api_layouts/`。
