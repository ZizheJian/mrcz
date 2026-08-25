import csv,os,re,statistics
from third_party.run_command.run_command import run_command

# XXX 比较正式双tree SZ3-3D-GLL与experimental五tree版本。
# XXX 五tree使用一棵角点tree和四棵预测层tree，因此只测试至少有四个预测层的element。
# XXX mixlay0和tc_channel0处理全部element；pb0和jicf0复用e18的四段100000-element样本，便于直接比较4tree。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
datasets=[
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"ranges":[(None,None)],"error_args":"-M REL 1e-4","max_cpus":1},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"ranges":[(None,None)],"error_args":"-M REL 1e-4","max_cpus":1},
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","ranges":[(162007,262007),(583513,683513),(949389,1049389),(1179079,1279079)],"error_args":"-M ABS 1.683629458334e-3","max_cpus":4},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","ranges":[(777593,877593),(1108559,1208559),(1245814,1345814),(2778116,2878116)],"error_args":"-M ABS 1.04476763947e-3","max_cpus":4},
]
methods=["SZ3-3D-GLL","SZ3-3D-GLL-experimental"]
field="UZ"
result_path=os.path.join(experiment_dir,"e19_sz3_five_huffman_tree_results.csv")

def parse_number(output,pattern,name):
    match=re.search(pattern,output)
    if match is None:
        raise RuntimeError(f"无法从输出中解析{name}")
    return float(match.group(1))

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup_method_outputs(data_path,method):
    remove_if_exists(f"{data_path}_{field}.{method}")
    remove_if_exists(f"{data_path}.{method}.bin")

rows=[]
fieldnames=["Dataset","Sample","Element Begin","Element End","Method","BR","PSNR","Max Error","Compression Throughput","Decompression Throughput"]
with open(result_path,"w",newline="") as result_file:
    writer=csv.DictWriter(result_file,fieldnames=fieldnames)
    writer.writeheader()
    for dataset in datasets:
        data_path=dataset["path"]
        auxiliary_paths=[f"{data_path}_header.mrcz",f"{data_path}_element_map.mrcz"]
        auxiliary_existed=[os.path.exists(path) for path in auxiliary_paths]
        for method in methods:
            compressed_path=f"{data_path}_{field}.{method}"
            method_bin_path=f"{data_path}.{method}.bin"
            if os.path.exists(compressed_path) or os.path.exists(method_bin_path):
                raise RuntimeError(f"拒绝覆盖已有实验文件: {compressed_path}或{method_bin_path}")
        try:
            for sample_id,(begin,end) in enumerate(dataset["ranges"],1):
                element_args="" if begin is None else f" -elements [{begin},{end})"
                for method in methods:
                    base_args=f"-format nekrs -i {data_path} {dataset['error_args']} -method {method} -max-cpus {dataset['max_cpus']}{element_args} -fields [{field}]"
                    compression_output=run_command(f"./install/bin/mrcz_compress {base_args}")
                    decompression_output=run_command(f"./install/bin/mrcz_decompress {base_args}")
                    row={
                        "Dataset":dataset["name"],"Sample":sample_id,"Element Begin":"" if begin is None else begin,"Element End":"" if end is None else end,"Method":method,
                        "BR":parse_number(decompression_output,r"BR= ([0-9.eE+-]+)","BR"),
                        "PSNR":parse_number(decompression_output,r"PSNR= ([0-9.eE+-]+)","PSNR"),
                        "Max Error":parse_number(decompression_output,r"max_err= ([0-9.eE+-]+)","max error"),
                        "Compression Throughput":parse_number(compression_output,r"compression bandwidth per cpu= ([0-9.eE+-]+)","compression throughput"),
                        "Decompression Throughput":parse_number(decompression_output,r"decompression bandwidth per cpu= ([0-9.eE+-]+)","decompression throughput"),
                    }
                    rows.append(row)
                    writer.writerow(row)
                    result_file.flush()
                    cleanup_method_outputs(data_path,method)
        finally:
            for method in methods:
                cleanup_method_outputs(data_path,method)
            for path,existed in zip(auxiliary_paths,auxiliary_existed):
                if not existed:
                    remove_if_exists(path)

for dataset in datasets:
    dataset_rows=[row for row in rows if row["Dataset"]==dataset["name"]]
    formal=[row for row in dataset_rows if row["Method"]=="SZ3-3D-GLL"]
    experimental=[row for row in dataset_rows if row["Method"]=="SZ3-3D-GLL-experimental"]
    br_change=(sum(row["BR"] for row in experimental)/sum(row["BR"] for row in formal)-1)*100
    formal_compression=statistics.harmonic_mean(row["Compression Throughput"] for row in formal)
    experimental_compression=statistics.harmonic_mean(row["Compression Throughput"] for row in experimental)
    formal_decompression=statistics.harmonic_mean(row["Decompression Throughput"] for row in formal)
    experimental_decompression=statistics.harmonic_mean(row["Decompression Throughput"] for row in experimental)
    print(f"{dataset['name']}: BR change={br_change:.4f}%, compression throughput={formal_compression:.4f}->{experimental_compression:.4f} MB/s/cpu, decompression throughput={formal_decompression:.4f}->{experimental_decompression:.4f} MB/s/cpu")

print(f"结果已写入{result_path}")
