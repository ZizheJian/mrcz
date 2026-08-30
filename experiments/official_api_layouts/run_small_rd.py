# XXX 比较所有官方API baseline以及MRCZ、SZ3强化版在小数据上的RD表现。
# XXX 每个方法使用一个覆盖完整field的partition和一个CPU，必要的逐element bitstream长度表计入BR。
# XXX 脚本支持断点续跑；已经存在于数据旁边的逐field CSV中的dataset/method/EB组合不会重复执行。

import csv,glob,os,re,subprocess,sys

project_root=os.path.abspath(os.path.join(os.path.dirname(__file__),"../.."))
sys.path.insert(0,project_root)
from third_party.csv_tools.py_codes.change_string_type import change_string_type
from third_party.csv_tools.py_codes.write_csv_item import write_csv_item
experiment_dir=os.path.dirname(os.path.abspath(__file__))
log_path=os.path.join(experiment_dir,"raw.log")
ansi_pattern=re.compile(r"\x1b\[[0-9;]*m")
methods=[
    "SZ3-1D-uniform","SZ3-2D-uniform","SZ3-3D-uniform-fused","SZ3-3D-uniform-per-element","SZ3-4D-uniform",
    "SZ2-1D-uniform","SZ2-2D-uniform","SZ2-3D-uniform-fused","SZ2-3D-uniform-per-element","SZ2-4D-uniform",
    "QoZ-1D-uniform","QoZ-2D-uniform","QoZ-3D-uniform-fused","QoZ-3D-uniform-per-element","QoZ-4D-uniform",
    "ZFP-1D-uniform","ZFP-2D-uniform","ZFP-3D-uniform-fused","ZFP-3D-uniform-per-element","ZFP-4D-uniform",
    "SPERR-2D-uniform","SPERR-3D-uniform-fused","SPERR-3D-uniform-per-element",
    "MRCZ-3D-GLL","SZ3-3D-GLL-optimized"]
result_method_names={method:method.replace("-uniform","") for method in methods if "-uniform" in method}
result_method_names.update({"MRCZ-3D-GLL":"MRCZ","SZ3-3D-GLL-optimized":"FDWZ"})
eb_list=[1e-2,1e-3,1e-4,1e-5,1e-6]
datasets=[
    {"name":"mixlay0","format":"nekrs","data":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"fields":["UX","UY","P","T","S00","S01"]},
    {"name":"test0","format":"nekrs","data":os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),"fields":["UZ","T","S00","S01"]},
    {"name":"tc_channel0","format":"nekrs","data":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"fields":["UX","UY","UZ"]},
    {"name":"mfem_ex0","format":"mfem","data":os.path.join(project_root,"build/mfem_build/examples/ex0examples/sol.gf"),"mesh":os.path.join(project_root,"build/mfem_build/examples/ex0examples/mesh.mesh"),"fields":["no_name_field"]},
    {"name":"mfem_ex1","format":"mfem","data":os.path.join(project_root,"build/mfem_build/examples/ex1examples/sol.gf"),"mesh":os.path.join(project_root,"build/mfem_build/examples/ex1examples/refined.mesh"),"fields":["no_name_field"]},
    {"name":"nektar_cylinder","format":"nektar++","data":os.path.join(project_root,"third_party/flow-stability-cylinder/tutorial/base/Cylinder-Base.vtu"),"fields":["u","v","p"]}]
# XXX 官方SZ2在20480次5x5x1逐element调用、REL=1e-4时于其内部分配器触发堆断言，原始失败保留在日志中。
known_failures={("mfem_ex1","SZ2-3D-uniform-per-element","0.0001")}

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def format_eb(value):
    return change_string_type(float(value))

def rate_distortion_path(dataset,field):
    return f"{dataset['data']}_{field}.csv"

def throughput_path(dataset,field):
    return f"{dataset['data']}_{field}_throughput.csv"

def result_method_name(method):
    return result_method_names.get(method,method)

def cleanup(dataset,method):
    data_path=dataset["data"]
    for path in glob.glob(f"{glob.escape(data_path)}_*.{glob.escape(method)}"):
        remove_if_exists(path)
    remove_if_exists(f"{data_path}.{method}.bin")
    remove_if_exists(f"{data_path}_header.mrcz")
    remove_if_exists(f"{data_path}_element_map.mrcz")

def run_command(args,log_file):
    command=" ".join(args)
    log_file.write(f"\n$ {command}\n")
    result=subprocess.run(args,cwd=project_root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace")
    log_file.write(result.stdout)
    log_file.flush()
    return result.returncode,ansi_pattern.sub("",result.stdout)

def make_args(task,dataset,method,eb):
    args=[os.path.join(project_root,"install/bin",f"mrcz_{task}"),"-format",dataset["format"],"-i",dataset["data"]]
    if "mesh" in dataset:
        args.extend(["-mesh",dataset["mesh"]])
    args.extend(["-M","REL",str(eb),"-method",method,"-max-cpus","1","-partition-num","1"])
    if "fields" in dataset:
        args.extend(["-fields",f"[{','.join(dataset['fields'])}]"])
    return args

def parse_throughput(text,task):
    pattern=rf"Field= ([^,\n]+), {task} bandwidth per cpu= ([0-9.eE+-]+) MB/s"
    return {field.strip():float(value) for field,value in re.findall(pattern,text)}

def parse_decompression(text):
    number=r"(?:[0-9.eE+-]+|[-+]?inf|nan)"
    pattern=rf"Field= ([^,\n]+), abs_eb= ({number}), rel_eb= ({number}), BR= ({number}), CR= ({number}), max_err= ({number}) \([^)]*\), MSE= ({number}), PSNR= ({number}) dB"
    return {field.strip():{"Abs EB":float(abs_eb),"BR":float(br),"Max Error":float(max_error),"PSNR":float(psnr)} for field,abs_eb,rel_eb,br,cr,max_error,mse,psnr in re.findall(pattern,text)}

def load_completed():
    completed=set(known_failures)
    for dataset in datasets:
        available={}
        for field in dataset["fields"]:
            path=rate_distortion_path(dataset,field)
            if not os.path.isfile(path):
                continue
            with open(path,newline="") as f:
                reader=csv.DictReader(f)
                eb_name="REL EB" if "REL EB" in (reader.fieldnames or []) else "REL_EB"
                for row in reader:
                    eb=format_eb(row.get(eb_name,""))
                    for method in methods:
                        csv_method=result_method_name(method)
                        if row.get(f"{csv_method}_BR","")!="" and row.get(f"{csv_method}_PSNR","")!="":
                            available.setdefault((method,eb),set()).add(field)
        for (method,eb),fields in available.items():
            if set(dataset["fields"]).issubset(fields):
                completed.add((dataset["name"],method,eb))
    return completed

def write_results(dataset,method,eb,compression_throughput,decompression_throughput,results):
    csv_method=result_method_name(method)
    for field,result in results.items():
        row={"REL EB":eb}
        write_csv_item(rate_distortion_path(dataset,field),row,f"{csv_method}_BR",result["BR"])
        write_csv_item(rate_distortion_path(dataset,field),row,f"{csv_method}_PSNR",result["PSNR"])
        throughput_row={"REL EB":eb,"Max CPUs":1}
        if field in compression_throughput:
            write_csv_item(throughput_path(dataset,field),throughput_row,f"{csv_method} Compression Throughput",compression_throughput[field])
        if field in decompression_throughput:
            write_csv_item(throughput_path(dataset,field),throughput_row,f"{csv_method} Decompression Throughput",decompression_throughput[field])

def main():
    completed=load_completed()
    with open(log_path,"a") as log_file:
        for dataset in datasets:
            if not os.path.isfile(dataset["data"]):
                print(f"跳过不存在的数据: {dataset['data']}",flush=True)
                continue
            for eb in eb_list:
                for method in methods:
                    eb_string=format_eb(eb)
                    key=(dataset["name"],method,eb_string)
                    if key in completed:
                        continue
                    print(f"{dataset['name']} REL={eb:g} {method}",flush=True)
                    cleanup(dataset,method)
                    compression_status,compression_output=run_command(make_args("compress",dataset,method,eb),log_file)
                    if compression_status!=0:
                        print("  压缩失败，详见raw.log",flush=True)
                        cleanup(dataset,method)
                        continue
                    decompression_status,decompression_output=run_command(make_args("decompress",dataset,method,eb),log_file)
                    compression_throughput=parse_throughput(compression_output,"compression")
                    decompression_throughput=parse_throughput(decompression_output,"decompression")
                    results=parse_decompression(decompression_output)
                    if decompression_status!=0 or not results:
                        print("  解压失败，详见raw.log",flush=True)
                    else:
                        print(f"  完成，共{len(results)}个field",flush=True)
                        write_results(dataset,method,eb_string,compression_throughput,decompression_throughput,results)
                    cleanup(dataset,method)

if __name__=="__main__":
    main()
