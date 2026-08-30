import glob,os
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command
from third_party.csv_tools.py_codes.write_csv_item import write_csv_item

redbegin="\033[91m"
redend="\033[0m"
silent_build=False
build_max_cpus=get_max_cpus(8)
build_before_test=os.environ.get("MRCZ_T1_BUILD_BEFORE_TEST","true").lower() in ["1","true","yes"] # 多个slurm实验并行时，应该先完成build和install，再把这个设为False，再开始运行

# method可选:
#   "MRCZ": 正式MRCZ，element内变换由编译器自动向量化，系数QB使用32x32分块转置
#   "MRCZ-ablation": 使用param中的MRCZ消融开关
#   "FDWZ": 使用最大有效Huffman tree数，并在每个partition采样选择维度顺序和linear/cubic预测
#   "FDWZ-ablation": 使用param中的FDWZ消融开关
#   "SZ3-1D-uniform"/"SZ3-2D-uniform"/"SZ3-3D-uniform-fused"/"SZ3-3D-uniform-per-element"/"SZ3-4D-uniform": 官方SZ3 API的五种tensorization
#   "QoZ-1D-uniform"/"QoZ-2D-uniform"/"QoZ-3D-uniform-fused"/"QoZ-3D-uniform-per-element"/"QoZ-4D-uniform": 官方QoZ 2.0 q4 API的五种tensorization
#   "SZ2-1D-uniform"/"SZ2-2D-uniform"/"SZ2-3D-uniform-fused"/"SZ2-3D-uniform-per-element"/"SZ2-4D-uniform": 官方SZ2 API的五种tensorization
#   "ZFP-1D-uniform"/"ZFP-2D-uniform"/"ZFP-3D-uniform-fused"/"ZFP-3D-uniform-per-element"/"ZFP-4D-uniform": 官方ZFP API的五种tensorization
#   "SPERR-2D-uniform"/"SPERR-3D-uniform-fused"/"SPERR-3D-uniform-per-element": 官方SPERR API支持的三种tensorization
# 正式方法、外部baseline及其tensorization见README.md第5.1节。
formal_method_list=["MRCZ","FDWZ","SZ3-4D-uniform","SZ2-3D-uniform-fused","QoZ-4D-uniform","ZFP-3D-uniform-fused","SPERR-2D-uniform"]
method_list=os.environ.get("MRCZ_T1_METHODS",",".join(formal_method_list)).split(",")
official_uniform_methods={
    "SZ3-1D-uniform","SZ3-2D-uniform","SZ3-3D-uniform-fused","SZ3-3D-uniform-per-element","SZ3-4D-uniform",
    "SZ2-1D-uniform","SZ2-2D-uniform","SZ2-3D-uniform-fused","SZ2-3D-uniform-per-element","SZ2-4D-uniform",
    "QoZ-1D-uniform","QoZ-2D-uniform","QoZ-3D-uniform-fused","QoZ-3D-uniform-per-element","QoZ-4D-uniform",
    "ZFP-1D-uniform","ZFP-2D-uniform","ZFP-3D-uniform-fused","ZFP-3D-uniform-per-element","ZFP-4D-uniform",
    "SPERR-2D-uniform","SPERR-3D-uniform-fused","SPERR-3D-uniform-per-element"}
eb_mode="REL"
write_to_csv=True
cleanup_after_each_test=True
project_root=os.path.dirname(os.path.abspath(__file__))

# format="nekrs"
# data_path=os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001") #1600*8*8*1
# eb_list=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]
# fields=["UX","UY","P","T","S00","S01"]

# format="nekrs"
# data_path=os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001") #512*6*6*6
# fields=["UZ","T","S00","S01"]
# eb_list=[1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7,3.162e-8,1e-8,3.162e-9,1e-9,3.162e-10,1e-10]

# format="nekrs"
# data_path=os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001") #512*8*8*8 (double)
# fields=["UX","UY","UZ"]
# eb_list=[1e-2,3.162e-3,1e-3,3.162e-4,1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7,3.162e-8,1e-8]

# format="nekrs"
# data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001" #1947154*8*8*8 #CR= 26.078980036305357
# eb_list=[1e-3]
# MRCZ的相对误差为1.04e-4时，T的CR=26.072790
# SZ3-3D-uniform的相对误差为2e-3时，T的CR=26.078980036305357
# ZFP-4D-uniform的相对误差为3.16e-1时，T的CR=27.974553
# eb_list=[1e-1,1e-2,1e-3,1e-4]

test_case=os.environ.get("MRCZ_T1_DATASET","jicf0")
format="nekrs"
if test_case=="pb0":
    data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001" #1947154*8*8*8
    fields=["UX","UY","UZ","P","T"]
    eb_list=[1e-2,1e-3,1e-4,1e-5,1e-6]
elif test_case=="jicf0":
    data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000" #3240264*10*10*10
    fields=["UX","UY","UZ","P","S00"]
    eb_list=[1e-4,1e-5,1e-6,1e-7,1e-8]
else:
    raise ValueError(f"未知的MRCZ_T1_DATASET: {test_case}")
max_cpus_list=[128]
if "MRCZ_T1_FIELDS" in os.environ:
    fields=[field for field in os.environ["MRCZ_T1_FIELDS"].split(",") if field]
if "MRCZ_T1_EBS" in os.environ:
    eb_list=[float(eb) for eb in os.environ["MRCZ_T1_EBS"].split(",") if eb]

# format="mfem"
# data_path="/home/x-zjian1/mrcz/build/mfem_build/examples/ex0examples/sol.gf"
# mesh_path="/home/x-zjian1/mrcz/build/mfem_build/examples/ex0examples/mesh.mesh"
# eb_list=[1e-2,3.162e-3,1e-3,3.162e-4,1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7]

# format="mfem"
# data_path="/home/x-zjian1/mrcz/build/mfem_build/examples/ex1examples/sol.gf"
# mesh_path="/home/x-zjian1/mrcz/build/mfem_build/examples/ex1examples/refined.mesh"
# eb_list=[1e-2,3.162e-3,1e-3,3.162e-4,1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7]

# format="nektar++"
# data_path="/home/x-zjian1/mrcz/c_codes/third_party/flow-stability-cylinder/tutorial/base/Cylinder-Base.vtu"
# eb_list=[1e-1,1e-2,1e-3,1e-4,1e-5]
# eb_list=[1e1,1e0,1e-1,1e-2,1e-3]

def remove_file_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup_test_outputs(data_path,method,selected_fields=None):
    if selected_fields is None:
        compressed_paths=glob.glob(f"{glob.escape(data_path)}_*.{glob.escape(method)}")
    else:
        compressed_paths=[f"{data_path}_{field}.{method}" for field in selected_fields]

    for path in compressed_paths:
        remove_file_if_exists(path)
    remove_file_if_exists(f"{data_path}_header.mrcz")
    remove_file_if_exists(f"{data_path}_element_map.mrcz")
    remove_file_if_exists(f"{data_path}.{method}.bin")

def remove_ansi_color(text):
    return text.replace("\033[91m","").replace("\033[0m","")

def parse_rate_distortion(output):
    results={}
    for line in output.splitlines():
        line=remove_ansi_color(line).strip()
        if not line.startswith("Field= "):
            continue
        segments=line.split(",")
        if len(segments)!=8:
            continue
        field=segments[0].split("=",1)[1].strip()
        br=float(segments[3].split("=",1)[1].strip())
        psnr=float(segments[7].split("=",1)[1].strip().split()[0])
        results[field]=(br,psnr)
    return results

def parse_field_throughput(output,task):
    throughput={}
    marker=f"{task} bandwidth per cpu= "
    for line in output.splitlines():
        line=remove_ansi_color(line).strip()
        if not line.startswith("Field= "):
            continue
        segments=line.split(",")
        if len(segments)!=2 or not segments[1].strip().startswith(marker):
            continue
        field=segments[0].split("=",1)[1].strip()
        throughput[field]=float(segments[1].strip()[len(marker):].split()[0])
    return throughput

def result_method_name(method):
    if method in official_uniform_methods:
        return method.replace("-uniform","")
    return method

######## build sperr ######## sperr的CMakeLists写的有问题，需要手动安装并链接
sperr_source_path="third_party/sperr"
sperr_build_path="build/sperr_build"
sperr_install_path="install/sperr_install"
if build_before_test:
    os.makedirs(sperr_build_path,exist_ok=True)
    if not os.path.exists(os.path.join(sperr_build_path,"CMakeCache.txt")):
        run_command(f"cmake -S \"{sperr_source_path}\" -B \"{sperr_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{sperr_install_path}\" \
                    -DBUILD_SHARED_LIBS=OFF",silent=silent_build)
    run_command(f"cmake --build {sperr_build_path} --parallel {build_max_cpus}",silent=silent_build)
    run_command(f"cmake --install {sperr_build_path}",silent=silent_build)
######## build mfem ########
mfem_source_path="third_party/mfem"
mfem_build_path="build/mfem_build"
mfem_install_path="install/mfem_install"
if build_before_test:
    os.makedirs(mfem_build_path,exist_ok=True)
    if not os.path.exists(os.path.join(mfem_build_path,"CMakeCache.txt")):
        run_command(f"cmake -S \"{mfem_source_path}\" -B \"{mfem_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{mfem_install_path}\"",
                    silent=silent_build)
    run_command(f"cmake --build {mfem_build_path} --parallel {build_max_cpus}",silent=silent_build)
    run_command(f"cmake --install {mfem_build_path}",silent=silent_build)
######## build main ########
main_source_path="c_codes"
main_build_path="build"
main_install_path="install"
if build_before_test:
    os.makedirs(main_build_path,exist_ok=True)
    run_command(f"cmake -S \"{main_source_path}\" -B \"{main_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{main_install_path}\"",
                silent=silent_build)
    run_command(f"cmake --build {main_build_path} --parallel {build_max_cpus}",silent=silent_build)
    run_command(f"cmake --install {main_build_path}",silent=silent_build)
for method in method_list:
    for eb in eb_list:
        for max_cpus in max_cpus_list:
            command_list=["./install/bin/mrcz_compress","./install/bin/mrcz_decompress"]
            if format in ["nekrs","nektar++"]:
                command_list=[f"{command} -format {format} -i {data_path}" for command in command_list]
            elif format in ["mfem"]:
                command_list=[f"{command} -format {format} -i {data_path} -mesh {mesh_path}" for command in command_list] # type: ignore
            command_list=[f"{command} -M {eb_mode} {eb} -method {method}" for command in command_list]
            command_list=[f"{command} -max-cpus {max_cpus}" for command in command_list]
            if "partition_num" in locals():
                command_list=[f"{command} -partition-num {partition_num}" for command in command_list] # type: ignore
            if "elements" in locals():
                command_list=[f"{command} -elements {elements}" for command in command_list] # type: ignore
            if "fields" in locals():
                command_list=[f"{command} -fields [{','.join(fields)}]" for command in command_list] # type: ignore
            compression_output=run_command(command_list[0])
            decompression_output=run_command(command_list[1])
            if write_to_csv:
                compression_throughput=parse_field_throughput(compression_output,"compression")
                decompression_throughput=parse_field_throughput(decompression_output,"decompression")
                rate_distortion=parse_rate_distortion(decompression_output)
                csv_method_name=result_method_name(method)
                for field,(br,psnr) in rate_distortion.items():
                    rate_distortion_path=f"{data_path}_{field}.csv"
                    write_csv_item(rate_distortion_path,{f"{eb_mode} EB":eb},f"{csv_method_name}_BR",br)
                    write_csv_item(rate_distortion_path,{f"{eb_mode} EB":eb},f"{csv_method_name}_PSNR",psnr)

                    throughput_path=f"{data_path}_{field}_throughput.csv"
                    throughput_row={f"{eb_mode} EB":eb,"Max CPUs":max_cpus}
                    if "partition_num" in locals():
                        throughput_row["Partition Num"]=partition_num # type: ignore
                    if field in compression_throughput:
                        write_csv_item(throughput_path,throughput_row,
                            f"{csv_method_name} Compression Throughput",compression_throughput[field])
                    if field in decompression_throughput:
                        write_csv_item(throughput_path,throughput_row,
                            f"{csv_method_name} Decompression Throughput",decompression_throughput[field])
            if cleanup_after_each_test:
                cleanup_test_outputs(data_path,method,fields if "fields" in locals() else None)
