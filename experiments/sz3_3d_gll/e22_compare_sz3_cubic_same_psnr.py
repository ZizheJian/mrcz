import csv,math,os,re
from third_party.run_command.run_command import run_command

# XXX 为mixlay0和tc_channel0补三点rate-distortion曲线，在共同PSNR区间比较线性与三次GLL插值。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
datasets=[
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001")},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001")},
]
methods=["SZ3-3D-GLL","SZ3-3D-GLL-cubic"]
rel_ebs=[3.162e-4,1e-4,3.162e-5]
field="UZ"
result_path=os.path.join(experiment_dir,"e22_sz3_cubic_same_psnr_results.csv")

def parse_number(output,pattern,name):
    match=re.search(pattern,output)
    if match is None:
        raise RuntimeError(f"无法从输出中解析{name}")
    return float(match.group(1))

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup(data_path,method):
    remove_if_exists(f"{data_path}_{field}.{method}")
    remove_if_exists(f"{data_path}.{method}.bin")

def interpolate_log_br(points,target_psnr):
    points=sorted(points,key=lambda point:point[1])
    for (br0,psnr0),(br1,psnr1) in zip(points,points[1:]):
        if psnr0<=target_psnr<=psnr1:
            ratio=(target_psnr-psnr0)/(psnr1-psnr0)
            return math.exp(math.log(br0)+ratio*(math.log(br1)-math.log(br0)))
    raise RuntimeError("目标PSNR不在插值区间内")

rows=[]
with open(result_path,"w",newline="") as result_file:
    writer=csv.DictWriter(result_file,fieldnames=["Dataset","Method","REL EB","BR","PSNR","Max Error"])
    writer.writeheader()
    for dataset in datasets:
        data_path=dataset["path"]
        auxiliary_paths=[f"{data_path}_header.mrcz",f"{data_path}_element_map.mrcz"]
        auxiliary_existed=[os.path.exists(path) for path in auxiliary_paths]
        try:
            for method in methods:
                for rel_eb in rel_ebs:
                    args=f"-format nekrs -i {data_path} -M REL {rel_eb} -method {method} -max-cpus 1 -fields [{field}]"
                    run_command(f"./install/bin/mrcz_compress {args}")
                    output=run_command(f"./install/bin/mrcz_decompress {args}")
                    row={"Dataset":dataset["name"],"Method":method,"REL EB":rel_eb,"BR":parse_number(output,r"BR= ([0-9.eE+-]+)","BR"),
                        "PSNR":parse_number(output,r"PSNR= ([0-9.eE+-]+)","PSNR"),"Max Error":parse_number(output,r"max_err= ([0-9.eE+-]+)","max error")}
                    rows.append(row)
                    writer.writerow(row)
                    result_file.flush()
                    cleanup(data_path,method)
        finally:
            for method in methods:
                cleanup(data_path,method)
            for path,existed in zip(auxiliary_paths,auxiliary_existed):
                if not existed:
                    remove_if_exists(path)

for dataset in datasets:
    dataset_rows=[row for row in rows if row["Dataset"]==dataset["name"]]
    method_points={method:[(row["BR"],row["PSNR"]) for row in dataset_rows if row["Method"]==method] for method in methods}
    lower=max(min(psnr for br,psnr in points) for points in method_points.values())
    upper=min(max(psnr for br,psnr in points) for points in method_points.values())
    target_psnr=(lower+upper)/2
    linear_br=interpolate_log_br(method_points["SZ3-3D-GLL"],target_psnr)
    cubic_br=interpolate_log_br(method_points["SZ3-3D-GLL-cubic"],target_psnr)
    print(f"{dataset['name']}: target PSNR={target_psnr:.6f} dB,linear BR={linear_br:.9f},cubic BR={cubic_br:.9f},change={(cubic_br/linear_br-1)*100:.4f}%")

print(f"结果已写入{result_path}")
