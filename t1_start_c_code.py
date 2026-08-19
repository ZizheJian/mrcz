import glob,os,subprocess
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command
from third_party.csv_tools.py_codes.write_csv_item import write_csv_item

redbegin="\033[91m"
redend="\033[0m"
silent_build=False
build_max_cpus=get_max_cpus(8)
run_max_cpus=128

# method可选:
#   "mrcz"
#   "sz3-0": 把数据展成1维的，然后用SZ3的1维压缩
#   "sz3-1": 对每个element执行sz3的3维压缩，在block内共用huffman tree
#   "sz3-hohs": 对每个element执行sz3的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "sz2-0": 对每个element执行sz2的3维压缩，在block内共用huffman tree
#   "sz2-1": 对每个element执行sz2的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "zfp": 对每个element执行zfp的3维压缩
#   "sperr": 对每个element执行sperr的3维压缩
method_list=["mrcz"]
eb_mode="REL"
write_to_csv=True
cleanup_after_each_test=True
project_root=os.path.dirname(os.path.abspath(__file__))

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
    # 兼容并清理统一命名前由各压缩器生成的旧文件。
    remove_file_if_exists(f"{data_path}_header.{method}")
    remove_file_if_exists(f"{data_path}_element_map.{method}")
    remove_file_if_exists(f"{data_path}.{method}.bin")
    remove_file_if_exists(f"{data_path}.mrcz.bin")

# format="nekrs"
# data_path=os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001") #1600*8*8*1
# fields=["P"]
# eb_list=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]
# fields=["UX","UY","UZ","P","T"]

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
# SZ3-HOHS的相对误差为7.9e-4时，T的CR=26.052934653717053
# SZ3-1的相对误差为2e-3时，T的CR=26.078980036305357
# SZ2-1的相对误差为4.2e-3时，T的CR=26.051779
# SPERR的相对误差为4.52e-2时，T的CR=26.084365
# ZFP的相对误差为3.16e-1时，T的CR=27.974553
# eb_list=[1e-1,1e-2,1e-3,1e-4]

format="nekrs"
data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000" #3240264*10*10*10
fields=["UX","UY","UZ","P","S00"]
# eb_list=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]
eb_list=[1e-4,1e-5,1e-6,1e-7,1e-8]

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

######## build sperr ######## sperr的CMakeLists写的有问题，需要手动安装并链接
sperr_source_path="third_party/sperr"
sperr_build_path="build/sperr_build"
sperr_install_path="install/sperr_install"
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
os.makedirs(main_build_path,exist_ok=True)
run_command(f"cmake -S \"{main_source_path}\" -B \"{main_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{main_install_path}\"",
            silent=silent_build)
run_command(f"cmake --build {main_build_path} --parallel {build_max_cpus}",silent=silent_build)
run_command(f"cmake --install {main_build_path}",silent=silent_build)
for method in method_list:
    for eb in eb_list:
        if format=="nekrs" or format=="nektar++":
            command_list=["./install/bin/mrcz_compress","./install/bin/mrcz_decompress"]
            if format in ["nekrs","nektar++"]:
                command_list=[f"{command} -format {format} -i {data_path}" for command in command_list]
            elif format in ["mfem"]:
                command_list=[f"{command} -format {format} -i {data_path} -mesh {mesh_path}" for command in command_list] # type: ignore
            command_list=[f"{command} -M {eb_mode} {eb} -method {method}" for command in command_list]
            command_list=[f"{command} -max-cpus {run_max_cpus}" for command in command_list]
            if "fields" in locals():
                command_list=[f"{command} -fields [{','.join(fields)}]" for command in command_list]
            run_command(command_list[0])
            ret=run_command(command_list[1])
        if write_to_csv:
            ret=ret.splitlines()
            for line in ret:
                if not line.startswith("\033[91mField="):
                    continue
                segments=line.strip().split(",")
                if len(segments)==2:
                    field=segments[0].strip().split(" ")[1]
                    if segments[1].startswith("compression"):
                        compression_bandwidth=segments[1].strip().split("= ")[1].split(" ")[0]
                    else:
                        decompression_bandwidth=segments[1].strip().split("= ")[1].split(" ")[0]
                elif len(segments)==8:
                    field=segments[0].strip().split(" ")[1]
                    br=float(segments[3].strip().split(" ")[1])
                    psnr=float(segments[7].strip().split(" ")[1])
                    csv_method_name=method.upper()
                    write_csv_item(f"{data_path}_{field}.csv",{f"{eb_mode}_EB":eb},f"{csv_method_name}_BR",br)
                    write_csv_item(f"{data_path}_{field}.csv",{f"{eb_mode}_EB":eb},f"{csv_method_name}_PSNR",psnr)
        if cleanup_after_each_test:
            cleanup_test_outputs(data_path,method,fields if "fields" in locals() else None)
# elif method=="sz3-0":
#     import numpy as np
#     from py_funcs.read_nekrs_header import read_nekrs_header
#     # 编译sz3
#     os.chdir("c_codes/third_party/sz3")
#     os.makedirs("build",exist_ok=True)
#     os.chdir("build")
#     run_command("cmake -DCMAKE_INSTALL_PREFIX=../install ..",silent=silent_build)
#     run_command("make",silent=silent_build)
#     run_command("make install",silent=silent_build)
#     os.chdir("../../../..")

#     def sz3_0_compress(data_path,data,field):
#         with open(f"{data_path}_{field}","wb") as f2:
#             f2.write(data.tobytes())
#         for eb in eb_list:
#             command=f"./c_codes/third_party/sz3/install/bin/sz3 -i {data_path}_{field} -z {data_path}_{field}.sz3 -o {data_path}_{field}.sz3.bin "
#             if header.np_dtype==np.float32:
#                 command+=f"-f "
#             else:
#                 command+=f"-d "
#             command+=f"-M {eb_mode} {eb} -1 {data.shape[0]} -a "
#             print(command)
#             result=run_command(command,silent=True)
#             lines=result.splitlines()
#             for line in lines:
#                 if line.startswith("compression time = "):
#                     compression_time=float(line[len("compression time = "):])
#                 if line.startswith("decompression time = "):
#                     decompression_time=float(line[len("decompression time = "):].strip().split()[0])
#                 if len(line.split(","))==3 and line.split(",")[2].strip().startswith("range="):
#                     data_range=float(line.split(",")[2].strip()[len("range="):])
#                 if line.startswith("compression ratio = "):
#                     cr=float(line[len("compression ratio = "):])
#                 if line.startswith("Max absolute error = "):
#                     max_err=float(line.split(",")[0][len("Max absolute error = "):])
#                 if line.startswith("PSNR = "):
#                     psnr=float(line.split(",")[0][len("PSNR = "):])
#             compression_bandwidth=(header.file_element_num*header.element_size*header.word_size)/(compression_time*1024*1024)
#             decompression_bandwidth=(header.file_element_num*header.element_size*header.word_size)/(decompression_time*1024*1024)
#             if eb_mode=="REL":
#                 abs_eb=eb*data_range
#             else:
#                 abs_eb=eb
#             br=(8*header.word_size)/cr
#             mse=(data_range**2)/(10**(psnr/10))
#             print(f"{redbegin}Field: {field}, compression bandwidth= {compression_bandwidth} MB/s , decompression bandwidth= {decompression_bandwidth} MB/s{redend}")
#             print(f"{redbegin}Field: {field}, br= {br}, cr= {cr}, max_err= {max_err} ({max_err/abs_eb*100}%%), mse= {mse}, psnr= {psnr} dB{redend}")
#             add_csv("sz3-0",field,abs_eb,br,psnr)
#         os.remove(f"{data_path}_{field}")
#         os.remove(f"{data_path}_{field}.sz3")
#         os.remove(f"{data_path}_{field}.sz3.bin")
        
#     f=open(data_path,"rb")
#     header=read_nekrs_header(f)
#     # 跳过element map
#     f.seek(header.file_element_num*4,1)
#     # 跳过coordinate
#     if "X" in header.variables:
#         f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
#     # velocity
#     if "U" in header.variables:
#         if any("U" in field for field in fields):
#             data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*3*header.element_size)
#             data=data.reshape((header.file_element_num,3,header.element_size))
#             data=data.transpose((1,0,2)).reshape((3,-1))
#             if "UX" in fields:
#                 sz3_0_compress(data_path,data[0],"UX")
#             if "UY" in fields:
#                 sz3_0_compress(data_path,data[1],"UY")
#             if "UZ" in fields:
#                 sz3_0_compress(data_path,data[2],"UZ")
#         else:
#             f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
#     if "P" in header.variables:
#         if "P" in fields:
#             data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*header.element_size)
#             sz3_0_compress(data_path,data,"P")
#         else:
#             f.seek(header.file_element_num*header.element_size*header.word_size,1)
#     if "T" in header.variables:
#         if "T" in fields:
#             data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*header.element_size)
#             sz3_0_compress(data_path,data,"T")
#         else:
#             f.seek(header.file_element_num*header.element_size*header.word_size,1)
#     for sid in range(header.passive_scalar_num):
#         passive_scalar_name=f"S{sid:02d}"
#         if passive_scalar_name in fields:
#             data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*header.element_size)
#             sz3_0_compress(data_path,data,passive_scalar_name)
#         else:
#             f.seek(header.file_element_num*header.element_size*header.word_size,1)
