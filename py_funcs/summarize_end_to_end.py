import csv,math,os,statistics

methods=["MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
configs={
    "pb0":{"rel_eb":1e-4,"cpus":64,"value_num":996942848,"fields":["UX","UY","UZ","P","T"],"rd_dir":"pb0","rd_prefix":"pb0.f00001","transfer_csv":"pb0_transfer_rel_eb.csv"},
    "jicf0":{"rel_eb":1e-6,"cpus":128,"value_num":3240264000,"fields":["UX","UY","UZ","P","S00"],"rd_dir":"jicf0","rd_prefix":"jicf0.f00000","transfer_csv":"jicf0_transfer_rel_eb.csv"}}
fieldnames=["Dataset","REL EB","Method","Rank","Allocated CPUs","Active Backend CPUs","Fields","Raw Bytes","Compressed Bytes","Compression Ratio","Compression Time","Transfer Time Trial 1","Transfer Time Trial 2","Transfer Time Trial 3","Transfer Time Mean","Transfer Time Median","Transfer Time Std Dev","Decompression Time","End-to-end Time Min","End-to-end Time Median","End-to-end Time Max"]

def read_rows(path):
    with open(path,newline="") as input_file:
        return list(csv.DictReader(input_file))

def selected_throughput_row(path,rel_eb,cpus):
    matches=[row for row in read_rows(path) if math.isclose(float(row["REL EB"]),rel_eb,rel_tol=1e-12) and int(row["Max CPUs"])==cpus]
    if len(matches)!=1:
        raise RuntimeError(f"{path}应当恰好有一个REL EB={rel_eb:g}、Max CPUs={cpus}的记录，实际为{len(matches)}个")
    return matches[0]

def summarize_dataset(project_root,result_dir,dataset,config):
    raw_mib=config["value_num"]*4/(1024*1024)
    backend_times={method:{"compression":0.0,"decompression":0.0} for method in methods}
    for field in config["fields"]:
        throughput_path=os.path.join(project_root,"results","RD",config["rd_dir"],f"{config['rd_prefix']}_{field}_throughput.csv")
        row=selected_throughput_row(throughput_path,config["rel_eb"],config["cpus"])
        for method in methods:
            backend_times[method]["compression"]+=raw_mib/(float(row[f"{method} Compression Throughput"])*config["cpus"])
            backend_times[method]["decompression"]+=raw_mib/(float(row[f"{method} Decompression Throughput"])*config["cpus"])
    transfer_path=os.path.join(result_dir,config["transfer_csv"])
    transfer_rows=read_rows(transfer_path)
    summaries=[]
    for method in methods:
        rows=sorted((row for row in transfer_rows if row["Method"]==method),key=lambda row:int(row["Trial"]))
        if len(rows)!=3 or [int(row["Trial"]) for row in rows]!=[1,2,3]:
            raise RuntimeError(f"{transfer_path}中的{method}必须恰好包含Trial 1、2、3")
        if any(int(row["Files"])!=len(config["fields"]) or not math.isclose(float(row["REL EB"]),config["rel_eb"],rel_tol=1e-12) for row in rows):
            raise RuntimeError(f"{transfer_path}中的{method}元数据不一致")
        compressed_bytes={int(row["Bytes Transferred"]) for row in rows}
        if len(compressed_bytes)!=1:
            raise RuntimeError(f"{transfer_path}中的{method}三次传输字节数不一致")
        compressed_bytes=compressed_bytes.pop()
        speeds=[float(row["Effective Bytes Per Second"]) for row in rows]
        if any(not math.isfinite(speed) or speed<=0 for speed in speeds):
            raise RuntimeError(f"{transfer_path}中的{method}传输速度无效")
        transfer_times=[compressed_bytes/speed for speed in speeds]
        compression_time=backend_times[method]["compression"]
        decompression_time=backend_times[method]["decompression"]
        end_to_end_times=[compression_time+transfer_time+decompression_time for transfer_time in transfer_times]
        raw_bytes=config["value_num"]*4*len(config["fields"])
        summaries.append({"Dataset":dataset,"REL EB":config["rel_eb"],"Method":method,"Allocated CPUs":config["cpus"],"Active Backend CPUs":config["cpus"],"Fields":len(config["fields"]),"Raw Bytes":raw_bytes,"Compressed Bytes":compressed_bytes,"Compression Ratio":raw_bytes/compressed_bytes,"Compression Time":compression_time,"Transfer Time Trial 1":transfer_times[0],"Transfer Time Trial 2":transfer_times[1],"Transfer Time Trial 3":transfer_times[2],"Transfer Time Mean":statistics.mean(transfer_times),"Transfer Time Median":statistics.median(transfer_times),"Transfer Time Std Dev":statistics.stdev(transfer_times),"Decompression Time":decompression_time,"End-to-end Time Min":min(end_to_end_times),"End-to-end Time Median":statistics.median(end_to_end_times),"End-to-end Time Max":max(end_to_end_times)})
    direct_rows=sorted((row for row in transfer_rows if row["Method"]=="Direct"),key=lambda row:int(row["Trial"]))
    if direct_rows:
        raw_bytes=config["value_num"]*4*len(config["fields"])
        if len(direct_rows)!=3 or [int(row["Trial"]) for row in direct_rows]!=[1,2,3]:raise RuntimeError(f"{transfer_path}中的Direct必须恰好包含Trial 1、2、3")
        if any(int(row["Files"])!=1 or int(row["Bytes Transferred"])!=raw_bytes or not math.isclose(float(row["REL EB"]),config["rel_eb"],rel_tol=1e-12) for row in direct_rows):raise RuntimeError(f"{transfer_path}中的Direct元数据不一致")
        speed_text=[row["Effective Bytes Per Second"].strip() for row in direct_rows]
        if any(speed_text) and not all(speed_text):raise RuntimeError(f"{transfer_path}中的Direct传输速度只填写了一部分")
        if all(speed_text):
            speeds=[float(value) for value in speed_text]
            if any(not math.isfinite(speed) or speed<=0 for speed in speeds):raise RuntimeError(f"{transfer_path}中的Direct传输速度无效")
            transfer_times=[raw_bytes/speed for speed in speeds]
            summaries.append({"Dataset":dataset,"REL EB":config["rel_eb"],"Method":"Direct","Allocated CPUs":0,"Active Backend CPUs":0,"Fields":len(config["fields"]),"Raw Bytes":raw_bytes,"Compressed Bytes":raw_bytes,"Compression Ratio":1.0,"Compression Time":0.0,"Transfer Time Trial 1":transfer_times[0],"Transfer Time Trial 2":transfer_times[1],"Transfer Time Trial 3":transfer_times[2],"Transfer Time Mean":statistics.mean(transfer_times),"Transfer Time Median":statistics.median(transfer_times),"Transfer Time Std Dev":statistics.stdev(transfer_times),"Decompression Time":0.0,"End-to-end Time Min":min(transfer_times),"End-to-end Time Median":statistics.median(transfer_times),"End-to-end Time Max":max(transfer_times)})
        else:
            summaries.append({"Dataset":dataset,"REL EB":config["rel_eb"],"Method":"Direct","Allocated CPUs":0,"Active Backend CPUs":0,"Fields":len(config["fields"]),"Raw Bytes":raw_bytes,"Compressed Bytes":raw_bytes,"Compression Ratio":1.0,"Compression Time":0.0,"Transfer Time Trial 1":0.0,"Transfer Time Trial 2":0.0,"Transfer Time Trial 3":0.0,"Transfer Time Mean":0.0,"Transfer Time Median":0.0,"Transfer Time Std Dev":0.0,"Decompression Time":0.0,"End-to-end Time Min":0.0,"End-to-end Time Median":0.0,"End-to-end Time Max":0.0})
    summaries.sort(key=lambda row:row["End-to-end Time Median"])
    for rank,row in enumerate(summaries,1):
        row["Rank"]=rank
    return summaries

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    summaries=[]
    for dataset,config in configs.items():
        summaries.extend(summarize_dataset(project_root,result_dir,dataset,config))
    output_path=os.path.join(result_dir,"end_to_end_summary.csv")
    with open(output_path,"w",newline="") as output_file:
        writer=csv.DictWriter(output_file,fieldnames=fieldnames,lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)
    print(output_path)

if __name__=="__main__":
    main()
