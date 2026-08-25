import csv,os,re,statistics
from third_party.run_command.run_command import run_command

# XXX 比较SZ3-3D-GLL在每个预测层内固定采用六种X/Y/Z维度顺序时的真实bitstream、PSNR和程序内部throughput。
# XXX 所有方法都使用正式双Huffman-tree布局，只改变维度先后，避免混入tree数量的影响。
# XXX 三个小数据处理全部element；两个大数据各测试四段由固定随机种子产生的16384-element样本，用于观察不同partition是否偏好不同顺序。

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
orders=["XYZ","XZY","YXZ","YZX","ZXY","ZYX"]
methods=[f"SZ3-3D-GLL-order-{order}" for order in orders]
field="UZ"
result_path=os.path.join(experiment_dir,"e20_sz3_dimension_order_results.csv")

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
fieldnames=["Dataset","Sample","Element Begin","Element End","Order","BR","PSNR","Max Error","Compression Throughput","Decompression Throughput"]
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
                for order,method in zip(orders,methods):
                    base_args=f"-format nekrs -i {data_path} -M REL 1e-4 -method {method} -max-cpus {dataset['max_cpus']}{element_args} -fields [{field}]"
                    compression_output=run_command(f"./install/bin/mrcz_compress {base_args}")
                    decompression_output=run_command(f"./install/bin/mrcz_decompress {base_args}")
                    row={
                        "Dataset":dataset["name"],"Sample":sample_id,"Element Begin":"" if begin is None else begin,"Element End":"" if end is None else end,"Order":order,
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
    print(dataset["name"])
    for sample_id in sorted(set(row["Sample"] for row in dataset_rows)):
        sample_rows=[row for row in dataset_rows if row["Sample"]==sample_id]
        best=min(sample_rows,key=lambda row:row["BR"])
        worst=max(sample_rows,key=lambda row:row["BR"])
        psnrs=[row["PSNR"] for row in sample_rows]
        print(f"  sample {sample_id}: best={best['Order']} BR={best['BR']:.9f},worst={worst['Order']} BR={worst['BR']:.9f},BR spread={(worst['BR']/best['BR']-1)*100:.4f}%,PSNR spread={max(psnrs)-min(psnrs):.6f} dB")
    order_means={order:statistics.mean(row["BR"] for row in dataset_rows if row["Order"]==order) for order in orders}
    print("  mean BR: "+", ".join(f"{order}={order_means[order]:.9f}" for order in orders))

print(f"结果已写入{result_path}")
