import csv,os,random,re,statistics
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command

# XXX 当前比较正式双tree SZ3-3D-GLL与experimental四tree版本。
# XXX 随机抽取连续片段以保留原文件中的element顺序；打乱element会改变Zstd所见的序列。
# XXX 默认由run_e11.sbatch申请4 CPU和约7.6GB内存，每个数据抽取4段100000个element。

datasets=[
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","element_num":1947154,"abs_eb":1.683629458334e-3},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","element_num":3240264,"abs_eb":1.04476763947e-3},
]
methods=["SZ3-3D-GLL","SZ3-3D-GLL-experimental"]
field="UZ"
sample_num=4
sample_element_num=100000
random_seed=20260821
max_cpus=get_max_cpus(min(4,int(os.environ.get("SLURM_CPUS_PER_TASK","4"))))
result_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"e18_sz3_four_huffman_tree_large_results.csv")

def parse_number(output,pattern,name):
    match=re.search(pattern,output)
    if match is None:
        raise RuntimeError(f"无法从输出中解析{name}")
    return float(match.group(1))

def random_nonoverlapping_ranges(element_num,rng):
    ranges=[]
    while len(ranges)<sample_num:
        begin=rng.randrange(element_num-sample_element_num+1)
        end=begin+sample_element_num
        if all(end<=old_begin or begin>=old_end for old_begin,old_end in ranges):
            ranges.append((begin,end))
    return sorted(ranges)

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup_method_outputs(data_path,method):
    remove_if_exists(f"{data_path}_{field}.{method}")
    remove_if_exists(f"{data_path}.{method}.bin")

rng=random.Random(random_seed)
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
            ranges=random_nonoverlapping_ranges(dataset["element_num"],rng)
            print(f"{dataset['name']} ranges={ranges}")
            for sample_id,(begin,end) in enumerate(ranges,1):
                for method in methods:
                    base_args=f"-format nekrs -i {data_path} -M ABS {dataset['abs_eb']} -method {method} -max-cpus {max_cpus} -elements [{begin},{end}) -fields [{field}]"
                    compression_output=run_command(f"./install/bin/mrcz_compress {base_args}")
                    decompression_output=run_command(f"./install/bin/mrcz_decompress {base_args}")
                    row={
                        "Dataset":dataset["name"],"Sample":sample_id,"Element Begin":begin,"Element End":end,"Method":method,
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
    regular=[row for row in dataset_rows if row["Method"]=="SZ3-3D-GLL"]
    experimental=[row for row in dataset_rows if row["Method"]=="SZ3-3D-GLL-experimental"]
    br_change=(sum(row["BR"] for row in experimental)/sum(row["BR"] for row in regular)-1)*100
    regular_compression=statistics.harmonic_mean(row["Compression Throughput"] for row in regular)
    experimental_compression=statistics.harmonic_mean(row["Compression Throughput"] for row in experimental)
    regular_decompression=statistics.harmonic_mean(row["Decompression Throughput"] for row in regular)
    experimental_decompression=statistics.harmonic_mean(row["Decompression Throughput"] for row in experimental)
    print(f"{dataset['name']}: BR change={br_change:.4f}%, compression throughput={regular_compression:.4f}->{experimental_compression:.4f} MB/s/cpu, decompression throughput={regular_decompression:.4f}->{experimental_decompression:.4f} MB/s/cpu")

print(f"结果已写入{result_path}")
