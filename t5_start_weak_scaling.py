# 功能：在固定目标PSNR和完整RD确定的计算精度下测试正式方法的单节点weak scaling。
# 最大CPU处理完整数据；较少CPU处理固定全局partition网格中的多个区间，并按总原始字节数/总内部时间汇总。

import argparse,csv,glob,math,os,re,subprocess

project_root=os.path.dirname(os.path.abspath(__file__))
result_dir=os.path.join(project_root,"results","weak_scaling")
os.makedirs(result_dir,exist_ok=True)
methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
datasets={
    "pb0":{
        "data_path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001",
        "element_num":1947154,"element_size":512,"word_size":4,"max_cpus":64,"target_psnr":100.0,
        "fields":["UX","UY","UZ","P","T"],
        "precision":{"MRCZ":"double","FDWZ":"double","SZ3-4D":"double","SZ2-3D-fused":"double","QoZ-4D":"double","ZFP-3D-fused":"float","SPERR-2D":"double"},
        "ranges":{"UX":16.885440784,"UY":16.399092633,"UZ":16.836294966,"P":607.992369070,"T":1.584238466},
        "rel_eb":{
            "MRCZ":{"UX":6.18820155760789e-5,"UY":6.16202504472058e-5,"UZ":6.00099141450851e-5,"P":7.39421707059804e-5,"T":5.41583332061048e-5},
            "FDWZ":{"UX":1.94634703306144e-5,"UY":1.94424491898358e-5,"UZ":1.93190023602854e-5,"P":1.9820971362278e-5,"T":1.82275530426865e-5},
            "SZ3-4D":{"UX":1.78717125645538e-5,"UY":1.78674453292424e-5,"UZ":1.78567168288597e-5,"P":1.81123761371486e-5,"T":1.76458549384875e-5},
            "SZ2-3D-fused":{"UX":1.74685258451192e-5,"UY":1.7467352183969e-5,"UZ":1.74838056279523e-5,"P":1.74198220186186e-5,"T":1.74263082041909e-5},
            "QoZ-4D":{"UX":1.76750843862657e-5,"UY":1.76722425447876e-5,"UZ":1.76632451073545e-5,"P":1.78249436263498e-5,"T":1.74483636259059e-5},
            "ZFP-3D-fused":{"UX":3.73475695316773e-4,"UY":3.63749723475868e-4,"UZ":3.66864798548482e-4,"P":5.3357989972005e-4,"T":4.27366932251908e-4},
            "SPERR-2D":{"UX":2.45324593100814e-5,"UY":2.4528918670125e-5,"UZ":2.45010467406206e-5,"P":2.52839633058442e-5,"T":2.45011386435203e-5}},
        "abs_eb":{
            "MRCZ":{"UX":1.04490510960445e-3,"UY":1.01051619515239e-3,"UZ":1.01034461543099e-3,"P":4.49562755417074e-2,"T":8.57997147195563e-5},
            "FDWZ":{"UX":3.28649275718731e-4,"UY":3.18838525276513e-4,"UZ":3.25260422186615e-4,"P":1.205099933582e-2,"T":2.88767906712792e-5},
            "SZ3-4D":{"UX":3.01771744217442e-4,"UY":2.93009891069309e-4,"UZ":3.00640951655018e-4,"P":1.10121864771119e-2,"T":2.79552421590079e-5},
            "SZ2-3D-fused":{"UX":2.94963758741534e-4,"UY":2.86448726518142e-4,"UZ":2.94362508680415e-4,"P":1.05911188578776e-2,"T":2.76074277774506e-5},
            "QoZ-4D":{"UX":2.98451590756493e-4,"UY":2.89808742524815e-4,"UZ":2.97383604684177e-4,"P":1.08374297039236e-2,"T":2.76423688249153e-5},
            "ZFP-3D-fused":{"UX":6.3063017375346e-3,"UY":5.9651654105089e-3,"UZ":6.17664396100441e-3,"P":3.24412507318927e-1,"T":6.77051133169889e-4},
            "SPERR-2D":{"UX":4.14241388966268e-4,"UY":4.02252009458703e-4,"UZ":4.12506849900842e-4,"P":1.53724567497992e-2,"T":3.88156462998639e-5}}},
    "jicf":{
        "data_path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000",
        "element_num":3240264,"element_size":1000,"word_size":4,"max_cpus":128,"target_psnr":150.0,
        "fields":["UX","UY","UZ","P","S00"],
        "precision":{"MRCZ":"double","FDWZ":"double","SZ3-4D":"double","SZ2-3D-fused":"double","QoZ-4D":"double","ZFP-3D-fused":"double","SPERR-2D":"double"},
        "ranges":{"UX":9.3756573,"UY":8.9088866,"UZ":10.4476766,"P":40.3769319,"S00":1007636.1266886},
        "rel_eb":{
            "MRCZ":{"UX":2.27515184292434e-7,"UY":2.20136296910672e-7,"UZ":2.35990976578629e-7,"P":2.81007467480041e-7,"S00":5.18502936519852e-7},
            "FDWZ":{"UX":6.62930994918182e-8,"UY":6.69972345296007e-8,"UZ":6.91081024368415e-8,"P":7.24307797910966e-8,"S00":2.04507139533609e-7},
            "SZ3-4D":{"UX":5.75850311310915e-8,"UY":5.76160618944713e-8,"UZ":5.77599888683379e-8,"P":6.05325430390195e-8,"S00":2.16284527407402e-7},
            "SZ2-3D-fused":{"UX":5.50394266359839e-8,"UY":5.502572361184e-8,"UZ":5.49602600013848e-8,"P":5.49663677415258e-8,"S00":2.30603013832487e-7},
            "QoZ-4D":{"UX":5.68791272219437e-8,"UY":5.69207051453843e-8,"UZ":5.69993837685563e-8,"P":5.93331583414349e-8,"S00":1.76077978948221e-7},
            "ZFP-3D-fused":{"UX":1.60827378995462e-6,"UY":1.46859878370393e-6,"UZ":1.83233452297797e-6,"P":2.16692489319766e-6,"S00":5.88750609484843e-6},
            "SPERR-2D":{"UX":8.13448629175989e-8,"UY":8.13998817079611e-8,"UZ":8.16984381575428e-8,"P":8.62716244980366e-8,"S00":3.08232944079437e-7}},
        "abs_eb":{
            "MRCZ":{"UX":2.1331043984722e-6,"UY":1.96116930572111e-6,"UZ":2.46555740381169e-6,"P":1.13462193778331e-5,"S00":5.22462290631529e-1},
            "FDWZ":{"UX":6.21541382190092e-7,"UY":5.96870764937817e-7,"UZ":7.22019104699792e-7,"P":2.924532663089e-6,"S00":2.06068781959811e-1},
            "SZ3-4D":{"UX":5.39897517494945e-7,"UY":5.13294961756426e-7,"UZ":6.03457684115994e-7,"P":2.44411836802031e-6,"S00":2.17936103459468e-1},
            "SZ2-3D-fused":{"UX":5.16030802127477e-7,"UY":4.90217931740825e-7,"UZ":5.74207022346383e-7,"P":2.21937328708994e-6,"S00":2.32363927660885e-1},
            "QoZ-4D":{"UX":5.33279204356045e-7,"UY":5.07100107332265e-7,"UZ":5.95511128013165e-7,"P":2.39569089376403e-6,"S00":1.77422532702542e-1},
            "ZFP-3D-fused":{"UX":1.50786238991867e-5,"UY":1.30835800249163e-5,"UZ":1.91436385190891e-5,"P":8.74937788450568e-5,"S00":5.93246383726859},
            "SPERR-2D":{"UX":7.62661557830886e-7,"UY":7.2518231538964e-7,"UZ":8.53558860595107e-7,"P":3.4833835072596e-6,"S00":3.10586649890028e-1}}}}
ansi_pattern=re.compile(r"\x1b\[[0-9;]*m")
fieldnames=["Dataset","Target PSNR","Method","Field","CPUs","Trial","Group Index","Global Partition Begin","Global Partition End","Element Begin","Element End","Element Num","Partition Num","REL EB from full RD","ABS EB","Compute Precision","Raw Bytes","Value Num","Compression Time","Decompression Time","BR","MSE","PSNR using local range","PSNR using full range","Max Error"]
integer_columns=["CPUs","Trial","Group Index","Global Partition Begin","Global Partition End","Element Begin","Element End","Element Num","Partition Num","Raw Bytes","Value Num"]
float_columns=["Target PSNR","REL EB from full RD","ABS EB","Compression Time","Decompression Time","BR","MSE","PSNR using local range","PSNR using full range","Max Error"]

def extract(text,pattern,name):
    match=re.search(pattern,text)
    if match is None:
        raise RuntimeError(f"没有从输出中找到{name}")
    return float(match.group(1))

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup(config,method,field):
    data_path=config["data_path"]
    for path in glob.glob(f"{glob.escape(data_path)}_{glob.escape(field)}.{glob.escape(method)}"):
        remove_if_exists(path)
    remove_if_exists(f"{data_path}.{method}.bin")
    remove_if_exists(f"{data_path}_header.mrcz")
    remove_if_exists(f"{data_path}_element_map.mrcz")

def selected_groups(max_partitions,cpus):
    group_num=max_partitions//cpus
    if group_num<=4:
        group_indices=[]
        repeat=4//group_num
        for group_index in range(group_num):
            group_indices.extend([group_index]*repeat)
        return group_indices
    return sorted({round(index*(group_num-1)/3) for index in range(4)})

def group_element_range(config,cpus,group_index):
    max_partitions=config["max_cpus"]
    global_partition_begin=group_index*cpus
    global_partition_end=global_partition_begin+cpus
    element_begin=global_partition_begin*config["element_num"]//max_partitions
    element_end=global_partition_end*config["element_num"]//max_partitions
    return global_partition_begin,global_partition_end,element_begin,element_end

def run_command(args,log_file):
    command=" ".join(args)
    print(command,flush=True)
    log_file.write(f"\n$ {command}\n")
    result=subprocess.run(args,cwd=project_root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace")
    print(result.stdout,end="",flush=True)
    log_file.write(result.stdout)
    log_file.flush()
    if result.returncode!=0:
        raise subprocess.CalledProcessError(result.returncode,args,output=result.stdout)
    return ansi_pattern.sub("",result.stdout)

def make_args(config,task,method,field,cpus,element_begin,element_end,abs_eb):
    args=[os.path.join(project_root,"install/bin",f"mrcz_{task}"),"-format","nekrs","-i",config["data_path"],"-M","ABS",f"{abs_eb:.17g}","-method",method,"-max-cpus",str(cpus),"-partition-num",str(cpus),"-elements",f"[{element_begin},{element_end})","-fields",f"[{field}]"]
    if task=="compress":
        args.extend(["-compute-precision",config["precision"][method],"-constant-compression","false"])
    return args

def run_case(config,dataset,method,field,cpus,group_index,trial,log_file):
    global_begin,global_end,element_begin,element_end=group_element_range(config,cpus,group_index)
    abs_eb=config["abs_eb"][method][field]
    cleanup(config,method,field)
    compression_output=run_command(make_args(config,"compress",method,field,cpus,element_begin,element_end,abs_eb),log_file)
    decompression_output=run_command(make_args(config,"decompress",method,field,cpus,element_begin,element_end,abs_eb),log_file)
    compression_cpus=int(extract(compression_output,r"Max CPUs= ([0-9]+)","compression actual CPUs"))
    decompression_cpus=int(extract(decompression_output,r"Max CPUs= ([0-9]+)","decompression actual CPUs"))
    if compression_cpus!=cpus or decompression_cpus!=cpus:
        raise RuntimeError(f"请求{cpus}个CPU，但压缩/解压实际使用{compression_cpus}/{decompression_cpus}个CPU")
    compute_precision_match=re.search(r"Compute precision= (float|double)",compression_output)
    if compute_precision_match is None or compute_precision_match.group(1)!=config["precision"][method]:
        raise RuntimeError(f"请求{config['precision'][method]}计算精度，但压缩输出不一致")
    value_num=(element_end-element_begin)*config["element_size"]
    mse=extract(decompression_output,r"MSE= ([0-9.eE+-]+)","MSE")
    row={"Dataset":dataset,"Target PSNR":config["target_psnr"],"Method":method,"Field":field,"CPUs":cpus,"Trial":trial,"Group Index":group_index,"Global Partition Begin":global_begin,"Global Partition End":global_end,"Element Begin":element_begin,"Element End":element_end,"Element Num":element_end-element_begin,"Partition Num":cpus,"REL EB from full RD":config["rel_eb"][method][field],"ABS EB":abs_eb,"Compute Precision":config["precision"][method],"Raw Bytes":value_num*config["word_size"],"Value Num":value_num}
    row["Compression Time"]=extract(compression_output,r"Compression time= ([0-9.eE+-]+) seconds","compression time")
    row["Decompression Time"]=extract(decompression_output,r"Decompression time= ([0-9.eE+-]+) seconds","decompression time")
    row["BR"]=extract(compression_output,r"BR= ([0-9.eE+-]+) bits/value","BR")
    row["MSE"]=mse
    row["PSNR using local range"]=extract(decompression_output,r"PSNR= ([0-9.eE+-]+) dB","PSNR")
    row["PSNR using full range"]=20*math.log10(config["ranges"][field])-10*math.log10(mse)
    row["Max Error"]=extract(decompression_output,r"max_err= ([0-9.eE+-]+)","max error")
    cleanup(config,method,field)
    return row

def normalize_row(row):
    result=dict(row)
    for column in integer_columns:
        result[column]=int(result[column])
    for column in float_columns:
        result[column]=float(result[column])
    result.setdefault("Compute Precision","")
    return result

def row_key(row):
    return row["Method"],row["Field"],row["CPUs"],row["Trial"],row["Element Begin"],row["Element End"]

def write_raw_rows(path,rows):
    with open(path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def write_summaries(rows,config,dataset):
    summary=[]
    quality=[]
    cpu_counts=[1,2,4,8,16,32,64]+([128] if config["max_cpus"]==128 else [])
    for method in methods:
        method_rows=[row for row in rows if row["Method"]==method]
        baseline_rows=[row for row in method_rows if row["CPUs"]==1]
        if not baseline_rows:
            continue
        baseline_compression=sum(row["Raw Bytes"] for row in baseline_rows)/sum(row["Compression Time"] for row in baseline_rows)/1024/1024
        baseline_decompression=sum(row["Raw Bytes"] for row in baseline_rows)/sum(row["Decompression Time"] for row in baseline_rows)/1024/1024
        for cpus in cpu_counts:
            group=[row for row in method_rows if row["CPUs"]==cpus]
            if not group:
                continue
            compression=sum(row["Raw Bytes"] for row in group)/sum(row["Compression Time"] for row in group)/1024/1024
            decompression=sum(row["Raw Bytes"] for row in group)/sum(row["Decompression Time"] for row in group)/1024/1024
            value_num=sum(row["Value Num"] for row in group)
            summary.append({"Dataset":dataset,"Target PSNR":config["target_psnr"],"Method":method,"CPUs":cpus,"Samples":len(group),"Sampled Element Num":sum(row["Element Num"] for row in group),"Compression Throughput MB/s":compression,"Decompression Throughput MB/s":decompression,"Compression Throughput Per CPU MB/s":compression/cpus,"Decompression Throughput Per CPU MB/s":decompression/cpus,"Compression Efficiency":compression/cpus/baseline_compression,"Decompression Efficiency":decompression/cpus/baseline_decompression,"Aggregate BR":sum(row["BR"]*row["Value Num"] for row in group)/value_num})
            for field in config["fields"]:
                field_group=[row for row in group if row["Field"]==field]
                if not field_group:
                    continue
                field_values=sum(row["Value Num"] for row in field_group)
                mse=sum(row["MSE"]*row["Value Num"] for row in field_group)/field_values
                quality.append({"Dataset":dataset,"Target PSNR":config["target_psnr"],"Method":method,"Field":field,"CPUs":cpus,"REL EB from full RD":field_group[0]["REL EB from full RD"],"ABS EB":field_group[0]["ABS EB"],"Aggregate BR":sum(row["BR"]*row["Value Num"] for row in field_group)/field_values,"Aggregate MSE":mse,"PSNR using full range":20*math.log10(config["ranges"][field])-10*math.log10(mse),"Max Error":max(row["Max Error"] for row in field_group)})
    summary_path=os.path.join(result_dir,f"weak_{dataset}_summary.csv")
    quality_path=os.path.join(result_dir,f"weak_{dataset}_quality.csv")
    if summary:
        with open(summary_path,"w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(summary[0]))
            writer.writeheader()
            writer.writerows(summary)
    if quality:
        with open(quality_path,"w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(quality[0]))
            writer.writeheader()
            writer.writerows(quality)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",choices=sorted(datasets))
    parser.add_argument("--dry-run",action="store_true")
    parser.add_argument("--methods",default=",".join(methods))
    parser.add_argument("--append-existing",action="store_true")
    args=parser.parse_args()
    config=datasets[args.dataset]
    selected_methods=[method for method in args.methods.split(",") if method]
    unknown_methods=[method for method in selected_methods if method not in methods]
    if unknown_methods:
        raise ValueError(f"未知方法: {unknown_methods}")
    cpu_counts=[1,2,4,8,16,32,64]+([128] if config["max_cpus"]==128 else [])
    if args.dry_run:
        for cpus in cpu_counts:
            print(cpus,selected_groups(config["max_cpus"],cpus),[group_element_range(config,cpus,index) for index in selected_groups(config["max_cpus"],cpus)])
        return
    raw_path=os.path.join(result_dir,f"weak_{args.dataset}_raw.csv")
    log_path=os.path.join(result_dir,f"weak_{args.dataset}_raw.log")
    rows=[]
    if args.append_existing and os.path.isfile(raw_path):
        with open(raw_path,newline="") as f:
            rows=[normalize_row(row) for row in csv.DictReader(f)]
    write_raw_rows(raw_path,rows)
    with open(log_path,"a" if args.append_existing else "w") as log_file:
        existing_keys={row_key(row) for row in rows}
        for cpus in cpu_counts:
            groups=selected_groups(config["max_cpus"],cpus)
            for trial,group_index in enumerate(groups,1):
                for field_index,field in enumerate(config["fields"]):
                    ordered_methods=selected_methods if (trial+field_index+cpus)%2 else list(reversed(selected_methods))
                    for method in ordered_methods:
                        if method=="SZ2-3D-fused" and cpus!=1:
                            continue
                        _,_,element_begin,element_end=group_element_range(config,cpus,group_index)
                        key=method,field,cpus,trial,element_begin,element_end
                        if key in existing_keys:
                            continue
                        row=run_case(config,args.dataset,method,field,cpus,group_index,trial,log_file)
                        rows.append(row)
                        existing_keys.add(key)
                        write_raw_rows(raw_path,rows)
                        write_summaries(rows,config,args.dataset)

if __name__=="__main__":
    main()
