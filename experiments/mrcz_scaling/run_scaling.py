# XXX 测试当前MRCZ与SZ3强化版的strong scaling和weak scaling。
# XXX 所有吞吐量均来自压缩器内部计时；脚本运行时间不参与结果计算。
# XXX strong scaling固定element范围和partition布局，weak scaling固定每核element数和每个partition的element数。

import csv,glob,os,re,statistics,subprocess

project_root=os.path.abspath(os.path.join(os.path.dirname(__file__),"../.."))
experiment_dir=os.path.dirname(os.path.abspath(__file__))
data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
field="S00"
strong_eb_mode="REL"
strong_eb="1e-6"
weak_eb_mode="ABS"
weak_eb="5e-5"
methods=["MRCZ-3D-GLL","SZ3-3D-GLL-optimized"]
cpu_counts=[1,2,4,8,16,32,64,128]
repeat_num=3
strong_start=1350000
strong_element_num=512000
strong_partition_num=128
weak_start=300000
weak_element_num_per_cpu=20000
raw_path=os.path.join(experiment_dir,"raw_results.csv")
summary_path=os.path.join(experiment_dir,"summary.csv")
log_path=os.path.join(experiment_dir,"raw.log")
ansi_pattern=re.compile(r"\x1b\[[0-9;]*m")

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup(method):
    for path in glob.glob(f"{glob.escape(data_path)}_{glob.escape(field)}.{glob.escape(method)}"):
        remove_if_exists(path)
    remove_if_exists(f"{data_path}.{method}.bin")
    remove_if_exists(f"{data_path}_header.mrcz")
    remove_if_exists(f"{data_path}_element_map.mrcz")

def extract(text,pattern,name):
    match=re.search(pattern,text)
    if match is None:
        raise RuntimeError(f"没有从输出中找到{name}")
    return float(match.group(1))

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

def make_args(task,method,cpus,partition_num,start,end,eb_mode,eb):
    return [os.path.join(project_root,"install/bin",f"mrcz_{task}"),"-format","nekrs","-i",data_path,"-M",eb_mode,eb,"-method",method,"-max-cpus",str(cpus),"-partition-num",str(partition_num),"-elements",f"[{start},{end})","-fields",f"[{field}]"]

def run_case(mode,method,trial,cpus,log_file):
    if mode=="strong":
        start=strong_start
        element_num=strong_element_num
        partition_num=strong_partition_num
        eb_mode=strong_eb_mode
        eb=strong_eb
    else:
        start=weak_start
        element_num=weak_element_num_per_cpu*cpus
        partition_num=cpus
        eb_mode=weak_eb_mode
        eb=weak_eb
    end=start+element_num
    cleanup(method)
    compression_output=run_command(make_args("compress",method,cpus,partition_num,start,end,eb_mode,eb),log_file)
    decompression_output=run_command(make_args("decompress",method,cpus,partition_num,start,end,eb_mode,eb),log_file)
    row={"Mode":mode,"Method":method,"Trial":trial,"CPUs":cpus,"Partition Num":partition_num,"Start Element":start,"End Element":end,"Element Num":element_num,"EB Mode":eb_mode,"EB":eb}
    row["Compression Actual CPUs"]=int(extract(compression_output,r"Max CPUs= ([0-9]+)","compression actual CPUs"))
    row["Decompression Actual CPUs"]=int(extract(decompression_output,r"Max CPUs= ([0-9]+)","decompression actual CPUs"))
    if row["Compression Actual CPUs"]!=cpus or row["Decompression Actual CPUs"]!=cpus:
        raise RuntimeError(f"请求{cpus}个CPU，但压缩/解压实际使用{row['Compression Actual CPUs']}/{row['Decompression Actual CPUs']}个CPU")
    row["Compression Time"]=extract(compression_output,r"Compression time= ([0-9.eE+-]+) seconds","compression time")
    row["Decompression Time"]=extract(decompression_output,r"Decompression time= ([0-9.eE+-]+) seconds","decompression time")
    row["Compression Throughput Per CPU"]=extract(compression_output,r"compression bandwidth per cpu= ([0-9.eE+-]+) MB/s","compression throughput")
    row["Decompression Throughput Per CPU"]=extract(decompression_output,r"decompression bandwidth per cpu= ([0-9.eE+-]+) MB/s","decompression throughput")
    row["BR"]=extract(compression_output,r"BR= ([0-9.eE+-]+) bits/value","BR")
    row["PSNR"]=extract(decompression_output,r"PSNR= ([0-9.eE+-]+) dB","PSNR")
    row["Max Error"]=extract(decompression_output,r"max_err= ([0-9.eE+-]+)","max error")
    cleanup(method)
    return row

def write_summary(rows):
    summary=[]
    for mode in ["strong","weak"]:
        for method in methods:
            method_rows=[row for row in rows if row["Mode"]==mode and row["Method"]==method]
            if not method_rows:
                continue
            baseline_cpus=min(row["CPUs"] for row in method_rows)
            baseline_rows=[row for row in method_rows if row["CPUs"]==baseline_cpus]
            baseline_compression=statistics.median(row["Compression Throughput Per CPU"] for row in baseline_rows)
            baseline_decompression=statistics.median(row["Decompression Throughput Per CPU"] for row in baseline_rows)
            for cpus in cpu_counts:
                group=[row for row in method_rows if row["CPUs"]==cpus]
                if not group:
                    continue
                compression=statistics.median(row["Compression Throughput Per CPU"] for row in group)
                decompression=statistics.median(row["Decompression Throughput Per CPU"] for row in group)
                compression_speedup=cpus*compression/(baseline_cpus*baseline_compression)
                decompression_speedup=cpus*decompression/(baseline_cpus*baseline_decompression)
                summary.append({"Mode":mode,"Method":method,"CPUs":cpus,"Compression Throughput Per CPU Median":compression,"Decompression Throughput Per CPU Median":decompression,"Aggregate Compression Throughput":cpus*compression,"Aggregate Decompression Throughput":cpus*decompression,"Compression Speedup":compression_speedup,"Decompression Speedup":decompression_speedup,"Compression Efficiency":compression/baseline_compression,"Decompression Efficiency":decompression/baseline_decompression,"BR Median":statistics.median(row["BR"] for row in group),"PSNR Median":statistics.median(row["PSNR"] for row in group),"Max Error Maximum":max(row["Max Error"] for row in group)})
    if not summary:
        return
    with open(summary_path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

def main():
    rows=[]
    fieldnames=["Mode","Method","Trial","CPUs","Compression Actual CPUs","Decompression Actual CPUs","Partition Num","Start Element","End Element","Element Num","EB Mode","EB","Compression Time","Decompression Time","Compression Throughput Per CPU","Decompression Throughput Per CPU","BR","PSNR","Max Error"]
    with open(raw_path,"w",newline="") as raw_file,open(log_path,"w") as log_file:
        writer=csv.DictWriter(raw_file,fieldnames=fieldnames)
        writer.writeheader()
        for mode in ["strong","weak"]:
            for trial in range(1,repeat_num+1):
                ordered_cpus=cpu_counts if trial%2==1 else list(reversed(cpu_counts))
                for cpus in ordered_cpus:
                    ordered_methods=methods if (trial+cpus)%2==1 else list(reversed(methods))
                    for method in ordered_methods:
                        row=run_case(mode,method,trial,cpus,log_file)
                        writer.writerow(row)
                        raw_file.flush()
                        rows.append(row)
                        write_summary(rows)
    write_summary(rows)

if __name__=="__main__":
    main()
