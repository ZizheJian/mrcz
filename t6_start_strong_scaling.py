# 功能：在固定完整数据范围、ABS EB、计算精度和partition网格下测试单节点strong scaling。
# 所有计时均读取压缩器与解压器输出的backend内部时间，不使用进程wall time。

import argparse,csv,math,os,re,statistics
from t5_start_weak_scaling import cleanup,datasets,extract,methods,project_root,run_command

result_dir=os.path.join(project_root,"results","strong_scaling")
os.makedirs(result_dir,exist_ok=True)
parallel_methods=[method for method in methods if method!="SZ2-3D-fused"]
fieldnames=["Dataset","Target PSNR","Method","Field","CPUs","Trial","Element Begin","Element End","Element Num","Partition Num","REL EB from full RD","ABS EB","Compute Precision","Raw Bytes","Value Num","Compression Time","Decompression Time","BR","MSE","PSNR using full range","Max Error"]
integer_columns=["CPUs","Trial","Element Begin","Element End","Element Num","Partition Num","Raw Bytes","Value Num"]
float_columns=["Target PSNR","REL EB from full RD","ABS EB","Compression Time","Decompression Time","BR","MSE","PSNR using full range","Max Error"]

def make_args(config,task,method,field,cpus,partition_num,abs_eb,measure_phase_time):
    args=[os.path.join(project_root,"install/bin",f"mrcz_{task}"),"-format","nekrs","-i",config["data_path"],"-M","ABS",f"{abs_eb:.17g}","-method",method,"-max-cpus",str(cpus),"-partition-num",str(partition_num),"-elements",f"[0,{config['element_num']})","-fields",f"[{field}]"]
    if measure_phase_time:
        args.extend(["-measure-phase-time","true"])
    if task=="compress":
        args.extend(["-compute-precision",config["precision"][method],"-constant-compression","false"])
    return args

def extract_runtime_configuration(output,name):
    match=re.search(r"Max CPUs= ([0-9]+), partition num= ([0-9]+)",output)
    if match is None:
        raise RuntimeError(f"没有从{name}输出中找到CPU和partition数量")
    return int(match.group(1)),int(match.group(2))

def run_case(config,dataset,method,field,cpus,partition_num,trial,measure_phase_time,log_file):
    abs_eb=config["abs_eb"][method][field]
    cleanup(config,method,field)
    compression_output=run_command(make_args(config,"compress",method,field,cpus,partition_num,abs_eb,measure_phase_time),log_file)
    decompression_output=run_command(make_args(config,"decompress",method,field,cpus,partition_num,abs_eb,measure_phase_time),log_file)
    compression_cpus,compression_partitions=extract_runtime_configuration(compression_output,"compression")
    decompression_cpus,decompression_partitions=extract_runtime_configuration(decompression_output,"decompression")
    if compression_cpus!=cpus or decompression_cpus!=cpus:
        raise RuntimeError(f"请求{cpus}个CPU，但压缩/解压实际使用{compression_cpus}/{decompression_cpus}个CPU")
    if compression_partitions!=partition_num or decompression_partitions!=partition_num:
        raise RuntimeError(f"请求{partition_num}个partition，但压缩/解压实际使用{compression_partitions}/{decompression_partitions}个partition")
    compute_precision_match=re.search(r"Compute precision= (float|double)",compression_output)
    if compute_precision_match is None or compute_precision_match.group(1)!=config["precision"][method]:
        raise RuntimeError(f"请求{config['precision'][method]}计算精度，但压缩输出不一致")
    value_num=config["element_num"]*config["element_size"]
    mse=extract(decompression_output,r"MSE= ([0-9.eE+-]+)","MSE")
    max_error=extract(decompression_output,r"max_err= ([0-9.eE+-]+)","max error")
    if max_error>abs_eb*1.0001:
        raise RuntimeError(f"{method}/{field}/{cpus} CPU超过100.01%误差容差: {max_error/abs_eb:.9f}")
    row={"Dataset":dataset,"Target PSNR":config["target_psnr"],"Method":method,"Field":field,"CPUs":cpus,"Trial":trial,"Element Begin":0,"Element End":config["element_num"],"Element Num":config["element_num"],"Partition Num":partition_num,"REL EB from full RD":config["rel_eb"][method][field],"ABS EB":abs_eb,"Compute Precision":config["precision"][method],"Raw Bytes":value_num*config["word_size"],"Value Num":value_num}
    row["Compression Time"]=extract(compression_output,r"Compression time= ([0-9.eE+-]+) seconds","compression time")
    row["Decompression Time"]=extract(decompression_output,r"Decompression time= ([0-9.eE+-]+) seconds","decompression time")
    row["BR"]=extract(compression_output,r"BR= ([0-9.eE+-]+) bits/value","BR")
    row["MSE"]=mse
    row["PSNR using full range"]=20*math.log10(config["ranges"][field])-10*math.log10(mse)
    row["Max Error"]=max_error
    cleanup(config,method,field)
    return row

def normalize_strong_row(row):
    result=dict(row)
    for column in integer_columns:
        result[column]=int(result[column])
    for column in float_columns:
        result[column]=float(result[column])
    return result

def row_key(row):
    return row["Method"],row["Field"],row["CPUs"],row["Trial"]

def write_rows(path,columns,rows):
    with open(path,"w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=columns,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def aggregate_throughput(rows,column):
    return sum(row["Raw Bytes"] for row in rows)/sum(row[column] for row in rows)/1024/1024

def write_summaries(rows,config,dataset,result_tag):
    summary=[]
    quality=[]
    for method in methods:
        method_rows=[row for row in rows if row["Method"]==method]
        if not method_rows:
            continue
        baseline=[row for row in method_rows if row["CPUs"]==1]
        baseline_compression=aggregate_throughput(baseline,"Compression Time") if baseline else None
        baseline_decompression=aggregate_throughput(baseline,"Decompression Time") if baseline else None
        for cpus in sorted({row["CPUs"] for row in method_rows}):
            group=[row for row in method_rows if row["CPUs"]==cpus]
            compression=aggregate_throughput(group,"Compression Time")
            decompression=aggregate_throughput(group,"Decompression Time")
            trial_compression=[]
            trial_decompression=[]
            for trial in sorted({row["Trial"] for row in group}):
                trial_group=[row for row in group if row["Trial"]==trial]
                trial_compression.append(aggregate_throughput(trial_group,"Compression Time"))
                trial_decompression.append(aggregate_throughput(trial_group,"Decompression Time"))
            compression_speedup=compression/baseline_compression if baseline_compression is not None else ""
            decompression_speedup=decompression/baseline_decompression if baseline_decompression is not None else ""
            value_num=sum(row["Value Num"] for row in group)
            summary.append({"Dataset":dataset,"Target PSNR":config["target_psnr"],"Method":method,"CPUs":cpus,"Trials":len({row["Trial"] for row in group}),"Fields":len({row["Field"] for row in group}),"Partition Num":group[0]["Partition Num"],"Compression Throughput MB/s":compression,"Compression Throughput Std Dev MB/s":statistics.stdev(trial_compression) if len(trial_compression)>1 else 0,"Decompression Throughput MB/s":decompression,"Decompression Throughput Std Dev MB/s":statistics.stdev(trial_decompression) if len(trial_decompression)>1 else 0,"Compression Speedup":compression_speedup,"Decompression Speedup":decompression_speedup,"Compression Efficiency":compression_speedup/cpus if compression_speedup!="" else "","Decompression Efficiency":decompression_speedup/cpus if decompression_speedup!="" else "","Aggregate BR":sum(row["BR"]*row["Value Num"] for row in group)/value_num})
            for field in config["fields"]:
                field_group=[row for row in group if row["Field"]==field]
                if not field_group:
                    continue
                field_values=sum(row["Value Num"] for row in field_group)
                mse=sum(row["MSE"]*row["Value Num"] for row in field_group)/field_values
                quality.append({"Dataset":dataset,"Target PSNR":config["target_psnr"],"Method":method,"Field":field,"CPUs":cpus,"Trials":len(field_group),"Partition Num":field_group[0]["Partition Num"],"REL EB from full RD":field_group[0]["REL EB from full RD"],"ABS EB":field_group[0]["ABS EB"],"Aggregate BR":sum(row["BR"]*row["Value Num"] for row in field_group)/field_values,"Aggregate MSE":mse,"PSNR using full range":20*math.log10(config["ranges"][field])-10*math.log10(mse),"Max Error":max(row["Max Error"] for row in field_group)})
    tag=f"_{result_tag}" if result_tag else ""
    if summary:
        write_rows(os.path.join(result_dir,f"strong_{dataset}{tag}_summary.csv"),list(summary[0]),summary)
    if quality:
        write_rows(os.path.join(result_dir,f"strong_{dataset}{tag}_quality.csv"),list(quality[0]),quality)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",choices=sorted(datasets))
    parser.add_argument("--cpu-counts")
    parser.add_argument("--methods",default=",".join(parallel_methods))
    parser.add_argument("--fields")
    parser.add_argument("--partition-num",type=int)
    parser.add_argument("--repeats",type=int,default=1)
    parser.add_argument("--append-existing",action="store_true")
    parser.add_argument("--result-tag",default="")
    parser.add_argument("--measure-phase-time",action="store_true")
    parser.add_argument("--dry-run",action="store_true")
    args=parser.parse_args()
    if args.result_tag and re.fullmatch(r"[A-Za-z0-9_-]+",args.result_tag) is None:
        raise ValueError("--result-tag只能包含字母、数字、下划线和连字符")
    config=datasets[args.dataset]
    available_cpu_counts=[1,2,4,8,16,32,64]+([128] if config["max_cpus"]==128 else [])
    cpu_counts=available_cpu_counts if args.cpu_counts is None else [int(value) for value in args.cpu_counts.split(",") if value]
    selected_methods=[method for method in args.methods.split(",") if method]
    selected_fields=config["fields"] if args.fields is None else [field for field in args.fields.split(",") if field]
    partition_num=config["max_cpus"] if args.partition_num is None else args.partition_num
    if not cpu_counts or any(cpus not in available_cpu_counts for cpus in cpu_counts):
        raise ValueError(f"CPU数量必须来自{available_cpu_counts}")
    if not selected_methods or any(method not in methods for method in selected_methods):
        raise ValueError(f"方法必须来自{methods}")
    if "SZ2-3D-fused" in selected_methods and any(cpus>1 for cpus in cpu_counts):
        raise ValueError("SZ2官方实现使用进程级全局状态，不能执行多CPU strong scaling")
    if not selected_fields or any(field not in config["fields"] for field in selected_fields):
        raise ValueError(f"field必须来自{config['fields']}")
    if partition_num<=0 or partition_num>config["element_num"]:
        raise ValueError("--partition-num超出有效范围")
    if args.repeats<=0:
        raise ValueError("--repeats必须为正整数")
    case_count=len(cpu_counts)*len(selected_methods)*len(selected_fields)*args.repeats
    if args.dry_run:
        print(f"dataset={args.dataset}, CPU counts={cpu_counts}, methods={selected_methods}, fields={selected_fields}, partition num={partition_num}, repeats={args.repeats}, cases={case_count}")
        return
    tag=f"_{args.result_tag}" if args.result_tag else ""
    raw_path=os.path.join(result_dir,f"strong_{args.dataset}{tag}_raw.csv")
    log_path=os.path.join(result_dir,f"strong_{args.dataset}{tag}_raw.log")
    rows=[]
    if args.append_existing and os.path.isfile(raw_path):
        with open(raw_path,newline="") as file:
            rows=[normalize_strong_row(row) for row in csv.DictReader(file)]
    write_rows(raw_path,fieldnames,rows)
    with open(log_path,"a" if args.append_existing else "w") as log_file:
        existing_keys={row_key(row) for row in rows}
        for cpus in cpu_counts:
            for trial in range(1,args.repeats+1):
                for field_index,field in enumerate(selected_fields):
                    ordered_methods=selected_methods if (cpus+trial+field_index)%2 else list(reversed(selected_methods))
                    for method in ordered_methods:
                        key=method,field,cpus,trial
                        if key in existing_keys:
                            print(f"跳过已有case: {key}",flush=True)
                            continue
                        row=run_case(config,args.dataset,method,field,cpus,partition_num,trial,args.measure_phase_time,log_file)
                        rows.append(row)
                        existing_keys.add(key)
                        write_rows(raw_path,fieldnames,rows)
                        print(f"完成{len(rows)}/{case_count}: {key}",flush=True)
    write_summaries(rows,config,args.dataset,args.result_tag)
    print(raw_path)

if __name__=="__main__":
    main()
