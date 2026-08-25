import csv,os,re
from third_party.run_command.run_command import run_command

# XXX 这是生成e16四tree历史结果时的脚本；当前experimental已改为五tree，新实验使用e19。
# XXX test0只有三个预测层，因此当前不能直接重跑本脚本。
# XXX 三个小数据处理全部element；pb0和jicf0各处理三个4096-element片段。默认单CPU，可直接在登录节点运行。
# XXX 此脚本此前也用于生成已保留的e12单tree和e13三tree历史结果；当前版本输出四tree的e16结果。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
datasets=(
    {"name":"test0","path":os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),"ranges":((None,None),),"error_args":"-M REL 1e-4"},
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"ranges":((None,None),),"error_args":"-M REL 1e-4"},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"ranges":((None,None),),"error_args":"-M REL 1e-4"},
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","ranges":((666618,670714),(955020,959116),(1398182,1402278)),"error_args":"-M ABS 1.683629458334e-3"},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","ranges":((294168,298264),(1139024,1143120),(2674406,2678502)),"error_args":"-M ABS 1.04476763947e-3"},
)
methods=("SZ3-3D-GLL","SZ3-3D-GLL-experimental")
field="UZ"
result_path=os.path.join(experiment_dir,"e16_sz3_four_huffman_tree_results.csv")

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
fieldnames=("Dataset","Sample","Element Begin","Element End","Method","BR","PSNR","Max Error")
with open(result_path,"w",newline="") as result_file:
    writer=csv.DictWriter(result_file,fieldnames=fieldnames)
    writer.writeheader()
    for dataset in datasets:
        data_path=dataset["path"]
        auxiliary_paths=(f"{data_path}_header.mrcz",f"{data_path}_element_map.mrcz")
        auxiliary_existed=tuple(os.path.exists(path) for path in auxiliary_paths)
        try:
            for sample_id,(begin,end) in enumerate(dataset["ranges"],1):
                element_args="" if begin is None else f" -elements [{begin},{end})"
                for method in methods:
                    base_args=f"-format nekrs -i {data_path} {dataset['error_args']} -method {method} -max-cpus 1{element_args} -fields [{field}]"
                    run_command(f"./install/bin/mrcz_compress {base_args}")
                    output=run_command(f"./install/bin/mrcz_decompress {base_args}")
                    row={"Dataset":dataset["name"],"Sample":sample_id,"Element Begin":"" if begin is None else begin,"Element End":"" if end is None else end,"Method":method,
                        "BR":parse_number(output,r"BR= ([0-9.eE+-]+)","BR"),"PSNR":parse_number(output,r"PSNR= ([0-9.eE+-]+)","PSNR"),
                        "Max Error":parse_number(output,r"max_err= ([0-9.eE+-]+)","max error")}
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
    print(f"{dataset['name']}: aggregate BR change={br_change:.4f}%")

print(f"结果已写入{result_path}")
