import csv,os,re,statistics
from third_party.run_command.run_command import run_command

# XXX 比较正式线性SZ3-3D-GLL和四点GLL坐标Lagrange三次插值版本。
# XXX 左右外点都存在时用四点GLL Lagrange三次插值；只存在一个外点时用它和两个近点做三点二次插值；没有外点时保持线性。
# XXX 三个小数据处理全部element；pb0和jicf0复用e20的四段16384-element固定随机样本。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
sample_size=16384
datasets=[
    {"name":"test0","path":os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),"ranges":[(None,None)],"max_cpus":1},
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"ranges":[(None,None)],"max_cpus":1},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"ranges":[(None,None)],"max_cpus":1},
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","ranges":[(1335789,1335789+sample_size),(1840664,1840664+sample_size),(1045503,1045503+sample_size),(1596147,1596147+sample_size)],"max_cpus":4},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","ranges":[(578353,578353+sample_size),(1229283,1229283+sample_size),(244351,244351+sample_size),(1273638,1273638+sample_size)],"max_cpus":4},
]
methods=["SZ3-3D-GLL","SZ3-3D-GLL-cubic"]
field="UZ"
result_path=os.path.join(experiment_dir,"e21_sz3_cubic_results.csv")

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
                    base_args=f"-format nekrs -i {data_path} -M REL 1e-4 -method {method} -max-cpus {dataset['max_cpus']}{element_args} -fields [{field}]"
                    compression_output=run_command(f"./install/bin/mrcz_compress {base_args}")
                    decompression_output=run_command(f"./install/bin/mrcz_decompress {base_args}")
                    row={
                        "Dataset":dataset["name"],"Sample":sample_id,"Element Begin":"" if begin is None else begin,"Element End":"" if end is None else end,"Method":method,
                        "BR":parse_number(decompression_output,r"BR= ([0-9.eE+-]+)","BR"),"PSNR":parse_number(decompression_output,r"PSNR= ([0-9.eE+-]+)","PSNR"),
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
    linear=[row for row in dataset_rows if row["Method"]=="SZ3-3D-GLL"]
    cubic=[row for row in dataset_rows if row["Method"]=="SZ3-3D-GLL-cubic"]
    br_change=(sum(row["BR"] for row in cubic)/sum(row["BR"] for row in linear)-1)*100
    psnr_changes=[cubic_row["PSNR"]-linear_row["PSNR"] for linear_row,cubic_row in zip(linear,cubic)]
    compression_change=(statistics.harmonic_mean(row["Compression Throughput"] for row in cubic)/statistics.harmonic_mean(row["Compression Throughput"] for row in linear)-1)*100
    decompression_change=(statistics.harmonic_mean(row["Decompression Throughput"] for row in cubic)/statistics.harmonic_mean(row["Decompression Throughput"] for row in linear)-1)*100
    print(f"{dataset['name']}: BR change={br_change:.4f}%,PSNR change=[{min(psnr_changes):.6f},{max(psnr_changes):.6f}] dB,compression throughput change={compression_change:.2f}%,decompression throughput change={decompression_change:.2f}%")
    for linear_row,cubic_row in zip(linear,cubic):
        print(f"  sample {linear_row['Sample']}: BR={linear_row['BR']:.9f}->{cubic_row['BR']:.9f} ({(cubic_row['BR']/linear_row['BR']-1)*100:.4f}%),PSNR={linear_row['PSNR']:.6f}->{cubic_row['PSNR']:.6f}")

print(f"结果已写入{result_path}")
