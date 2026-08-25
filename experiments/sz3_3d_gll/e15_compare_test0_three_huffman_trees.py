import csv,os,re
from third_party.run_command.run_command import run_command

# XXX 生成test0 UZ上SZ3-3D-uniform、正式双tree SZ3-3D-GLL和experimental四tree SZ3-3D-GLL的rate-distortion结果。
# XXX test0很小，默认单CPU并在登录节点运行；每次测试后删除bitstream和辅助文件。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
data_path=os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001")
field="UZ"
methods=("SZ3-3D-uniform","SZ3-3D-GLL","SZ3-3D-GLL-experimental")
eb_list=(1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7,3.162e-8,1e-8,3.162e-9,1e-9,3.162e-10,1e-10)
result_path=os.path.join(experiment_dir,"e17_test0_four_huffman_tree_rate_distortion.csv")

def parse_number(output,pattern,name):
    match=re.search(pattern,output)
    if match is None:
        raise RuntimeError(f"无法从输出中解析{name}")
    return float(match.group(1))

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup(method):
    remove_if_exists(f"{data_path}_{field}.{method}")
    remove_if_exists(f"{data_path}.{method}.bin")

auxiliary_paths=(f"{data_path}_header.mrcz",f"{data_path}_element_map.mrcz")
auxiliary_existed=tuple(os.path.exists(path) for path in auxiliary_paths)
fieldnames=("REL EB","Method","BR","PSNR","Max Error","Compression Throughput","Decompression Throughput")
with open(result_path,"w",newline="") as result_file:
    writer=csv.DictWriter(result_file,fieldnames=fieldnames)
    writer.writeheader()
    try:
        for rel_eb in eb_list:
            for method in methods:
                args=f"-format nekrs -i {data_path} -M REL {rel_eb} -method {method} -max-cpus 1 -fields [{field}]"
                compression_output=run_command(f"./install/bin/mrcz_compress {args}")
                decompression_output=run_command(f"./install/bin/mrcz_decompress {args}")
                writer.writerow({"REL EB":rel_eb,"Method":method,"BR":parse_number(decompression_output,r"BR= ([0-9.eE+-]+)","BR"),
                    "PSNR":parse_number(decompression_output,r"PSNR= ([0-9.eE+-]+)","PSNR"),"Max Error":parse_number(decompression_output,r"max_err= ([0-9.eE+-]+)","max error"),
                    "Compression Throughput":parse_number(compression_output,r"compression bandwidth per cpu= ([0-9.eE+-]+)","compression throughput"),
                    "Decompression Throughput":parse_number(decompression_output,r"decompression bandwidth per cpu= ([0-9.eE+-]+)","decompression throughput")})
                result_file.flush()
                cleanup(method)
    finally:
        for method in methods:
            cleanup(method)
        for path,existed in zip(auxiliary_paths,auxiliary_existed):
            if not existed:
                remove_if_exists(path)

print(f"结果已写入{result_path}")
