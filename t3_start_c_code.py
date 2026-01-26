import os,subprocess

redbegin="\033[91m"
redend="\033[0m"
silent_build=False

def run_command(command,silent=False):
    import pty,errno
    command_list=command.split()
    if silent:
        result=subprocess.run(command_list,capture_output=True,text=True,encoding="utf-8")
        if result.returncode!=0:
            print(result.stdout)
            print(result.stderr)
            exit()
        return result.stdout
    else:
        master,slave=pty.openpty()
        proc=subprocess.Popen(command_list,stdout=slave,stderr=slave,text=True,bufsize=1)
        os.close(slave)
        lines=[]
        with os.fdopen(master) as proc_stdout:
            try:
                for line in proc_stdout:
                    print(line,end="")
                    lines.append(line)
            except OSError as e:
                if e.errno!=errno.EIO:
                    raise e
        proc.wait()
        output="".join(lines)
        if proc.returncode!=0:
            print(output)
            exit()
        return output

# method可选:
#   "jzzz"
#   "sz3-0": 把数据展成1维的，然后用SZ3的1维压缩
#   "sz3-1": 对每个element执行sz3的3维压缩，在block内共用huffman tree
#   "sz3-2": 对每个element执行sz3的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "sz2-0": 对每个element执行sz2的3维压缩，在block内共用huffman tree
#   "sz2-1": 对每个element执行sz2的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "zfp": 对每个element执行zfp的3维压缩
#   "sperr": 对每个element执行sperr的3维压缩

method="jzzz"
data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
eb_mode="ABS"
eb=1e-6
fields=["P"]

if method=="jzzz" or method=="sz3-1" or method=="sz3-2" or method=="zfp" or method=="sz2-0" or method=="sz2-1" or method=="sperr":
    try:
        os.chdir("c_codes/third_party/sperr")
        os.makedirs("build",exist_ok=True)
        os.chdir("build")
        run_command("cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=OFF ..",silent=silent_build)
        run_command("make",silent=silent_build)
        run_command("make install",silent=silent_build)
        os.chdir("../../../..")
        os.makedirs("build",exist_ok=True)
        os.chdir("build")
        run_command("cmake -DCMAKE_INSTALL_PREFIX=../install ../c_codes",silent=silent_build)
        run_command("make",silent=silent_build)
        run_command("make install",silent=silent_build)
        os.chdir("..")
        os.system(f"./build/jzzz_compress -i {data_path} -M {eb_mode} {eb} -method {method} -fields [{','.join(fields)}]")
        os.system(f"./build/jzzz_decompress -i {data_path} -M {eb_mode} {eb} -method {method} -fields [{','.join(fields)}]")
        # os.system(f"./build/jzzz_compress -i {data_path} -M {eb_mode} {eb} -method {method}")
        # os.system(f"./build/jzzz_decompress -i {data_path} -M {eb_mode} {eb} -method {method}")
    except Exception as e:
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
        command=f"./c_codes/third_party/sz3/install/bin/sz3 -i {data_path}_{field} -z {data_path}_{field}.sz3 -o {data_path}_{field}.sz3.bin "
        command+=f"-f -M {eb_mode} {eb} -1 {data.shape[0]} -a "
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
        print(f"{redbegin}Field: {field}, br= {br} , cr= {cr} , max_err= {max_err} , mse= {mse}  ({max_err/abs_eb*100}%%) , psnr= {psnr}{redend}")
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
        if any("velocity" in field for field in fields):
            data=np.fromfile(f,dtype=header.np_dtype,count=header.file_element_num*3*header.element_size)
            data=data.reshape((header.file_element_num,3,header.element_size))
            data=data.transpose((1,0,2)).reshape((3,-1))
            if "velocityx" in fields:
                sz3_0_compress(data_path,data[0],"velocityx")
            if "velocityy" in fields:
                sz3_0_compress(data_path,data[1],"velocityy")
            if "velocityz" in fields:
                sz3_0_compress(data_path,data[2],"velocityz")
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