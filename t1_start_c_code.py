import os,subprocess
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command

redbegin="\033[91m"
redend="\033[0m"
silent_build=False
build_max_cpus=get_max_cpus(8)

# method可选:
#   "mrcz"
#   "sz3-0": 把数据展成1维的，然后用SZ3的1维压缩
#   "sz3-1": 对每个element执行sz3的3维压缩，在block内共用huffman tree
#   "sz3-hohs": 对每个element执行sz3的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "sz2-0": 对每个element执行sz2的3维压缩，在block内共用huffman tree
#   "sz2-1": 对每个element执行sz2的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "zfp": 对每个element执行zfp的3维压缩
#   "sperr": 对每个element执行sperr的3维压缩

method="mrcz"
eb_mode="REL"

format="nekrs"
data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/pySEMTools/examples/data/mixlay0.f00001" #1600*8*8*1
fields=["UX","UY","UZ","P","T"]
eb_list=[1e-1,3.162e-2,1e-2,3.162e-3,1e-3,3.162e-4,1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7]

# format="nekrs"
# data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/fld_pyutils/demos/data/test0.f00001" #512*6*6*6
# fields=["UZ","T","S00","S01"]
# eb_list=[1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7,3.162e-8,1e-8,3.162e-9,1e-9,3.162e-10,1e-10]

# format="nekrs"
# data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/pySEMTools/examples/data/tc_channel0.f00001" #512*8*8*8 (double)
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

# format="nekrs"
# data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000" #3240264*10*10*10
# fields=["UX","UY","UZ","P","S00"]
# eb_list=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]

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

if method in ["mrcz","sz3-hohs","sz3-1","sz2-1","sz2-0","zfp", "sperr"]:
    try:
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
        ######## build mrcz ########
        os.makedirs("build",exist_ok=True)
        os.chdir("build")
        run_command("cmake -DCMAKE_INSTALL_PREFIX=../install ../c_codes",silent=silent_build)
        run_command("make",silent=silent_build)
        run_command("make install",silent=silent_build)
        os.chdir("..")
        for eb in eb_list:
            if format=="nekrs" or format=="nektar++":
                if "fields" in locals():
                    run_command(f"./build/mrcz_compress -format {format} -i {data_path} -M {eb_mode} {eb} -method {method} -fields [{','.join(fields)}]")
                    # run_command(f"./build/mrcz_decompress -format {format} -i {data_path} -M {eb_mode} {eb} -method {method} -fields [{','.join(fields)}]")
                else:
                    os.system(f"./build/mrcz_compress -format {format} -i {data_path} -M {eb_mode} {eb} -method {method}")
                    ret=run_command(f"./build/mrcz_decompress -format {format} -i {data_path} -M {eb_mode} {eb} -method {method}")
            elif format=="mfem":
                if "fields" in locals():
                    os.system(f"./build/mrcz_compress -format {format} -i {data_path} -mesh {mesh_path} -M {eb_mode} {eb} -method {method} -fields [{','.join(fields)}]")
                    run_command(f"./build/mrcz_decompress -format {format} -i {data_path} -mesh {mesh_path} -M {eb_mode} {eb} -method {method} -fields [{','.join(fields)}]")
                else:
                    os.system(f"./build/mrcz_compress -format {format} -i {data_path} -mesh {mesh_path} -M {eb_mode} {eb} -method {method}")
                    run_command(f"./build/mrcz_decompress -format {format} -i {data_path} -mesh {mesh_path} -M {eb_mode} {eb} -method {method}")
            if "ret" in locals():
                ret=ret.splitlines()
                for line in ret:
                    segments=line.strip().split(",")
                    if len(segments)!=7:
                        continue
                    field=segments[0].strip().split(" ")[1]
                    abs_eb=float(segments[1].strip().split(" ")[1])
                    br=float(segments[2].strip().split(" ")[1])
                    psnr=float(segments[6].strip().split(" ")[1])
                    add_csv(method,field,abs_eb,br,psnr)
    except Exception as e:
        print(e)
        exit(0)
elif method=="sz3-0":
    import numpy as np
    from py_funcs.read_nekrs_header import read_nekrs_header
    # 编译sz3
    os.chdir("c_codes/third_party/sz3")
    os.makedirs("build",exist_ok=True)
    os.chdir("build")
    run_command("cmake -DCMAKE_INSTALL_PREFIX=../install ..",silent=silent_build)
    run_command("make",silent=silent_build)
    run_command("make install",silent=silent_build)
    os.chdir("../../../..")

    def sz3_0_compress(data_path,data,field):
        with open(f"{data_path}_{field}","wb") as f2:
            f2.write(data.tobytes())
        for eb in eb_list:
            command=f"./c_codes/third_party/sz3/install/bin/sz3 -i {data_path}_{field} -z {data_path}_{field}.sz3 -o {data_path}_{field}.sz3.bin "
            if header.np_dtype==np.float32:
                command+=f"-f "
            else:
                command+=f"-d "
            command+=f"-M {eb_mode} {eb} -1 {data.shape[0]} -a "
            print(command)
            result=run_command(command,silent=True)
            lines=result.splitlines()
            for line in lines:
                if line.startswith("compression time = "):
                    compression_time=float(line[len("compression time = "):])
                if line.startswith("decompression time = "):
                    decompression_time=float(line[len("decompression time = "):].strip().split()[0])
                if len(line.split(","))==3 and line.split(",")[2].strip().startswith("range="):
                    data_range=float(line.split(",")[2].strip()[len("range="):])
                if line.startswith("compression ratio = "):
                    cr=float(line[len("compression ratio = "):])
                if line.startswith("Max absolute error = "):
                    max_err=float(line.split(",")[0][len("Max absolute error = "):])
                if line.startswith("PSNR = "):
                    psnr=float(line.split(",")[0][len("PSNR = "):])
            compression_bandwidth=(header.file_element_num*header.element_size*header.word_size)/(compression_time*1024*1024)
            decompression_bandwidth=(header.file_element_num*header.element_size*header.word_size)/(decompression_time*1024*1024)
            if eb_mode=="REL":
                abs_eb=eb*data_range
            else:
                abs_eb=eb
            br=(8*header.word_size)/cr
            mse=(data_range**2)/(10**(psnr/10))
            print(f"{redbegin}Field: {field}, compression bandwidth= {compression_bandwidth} MB/s , decompression bandwidth= {decompression_bandwidth} MB/s{redend}")
            print(f"{redbegin}Field: {field}, br= {br}, cr= {cr}, max_err= {max_err} ({max_err/abs_eb*100}%%), mse= {mse}, psnr= {psnr} dB{redend}")
            add_csv("sz3-0",field,abs_eb,br,psnr)
        os.remove(f"{data_path}_{field}")
        os.remove(f"{data_path}_{field}.sz3")
        os.remove(f"{data_path}_{field}.sz3.bin")
        
    f=open(data_path,"rb")
    header=read_nekrs_header(f)
    # 跳过element map
    f.seek(header.file_element_num*4,1)
    # 跳过coordinate
    if "X" in header.variables:
        f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
    # velocity
    if "U" in header.variables:
        if any("U" in field for field in fields):
            data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*3*header.element_size)
            data=data.reshape((header.file_element_num,3,header.element_size))
            data=data.transpose((1,0,2)).reshape((3,-1))
            if "UX" in fields:
                sz3_0_compress(data_path,data[0],"UX")
            if "UY" in fields:
                sz3_0_compress(data_path,data[1],"UY")
            if "UZ" in fields:
                sz3_0_compress(data_path,data[2],"UZ")
        else:
            f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
    if "P" in header.variables:
        if "P" in fields:
            data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*header.element_size)
            sz3_0_compress(data_path,data,"P")
        else:
            f.seek(header.file_element_num*header.element_size*header.word_size,1)
    if "T" in header.variables:
        if "T" in fields:
            data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*header.element_size)
            sz3_0_compress(data_path,data,"T")
        else:
            f.seek(header.file_element_num*header.element_size*header.word_size,1)
    for sid in range(header.passive_scalar_num):
        passive_scalar_name=f"S{sid:02d}"
        if passive_scalar_name in fields:
            data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*header.element_size)
            sz3_0_compress(data_path,data,passive_scalar_name)
        else:
            f.seek(header.file_element_num*header.element_size*header.word_size,1)

