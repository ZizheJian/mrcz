import csv,math,os,statistics

project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
result_dir=os.path.join(project_root,"results","weak_scaling")
fields=["UX","UY","UZ","P","S00"]
methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
ranges={"UX":9.3756573,"UY":8.9088866,"UZ":10.4476766,"P":40.3769319,"S00":1007636.1266886}
expected_rows={1:640,2:320,4:160,8:80,16:40,32:20,64:20,128:20}
expected_repeats={1:1,2:1,4:1,8:1,16:1,32:1,64:2,128:4}
element_num=3240264
one_resume_identities={(126,127,"P"),(126,127,"S00")}|{(127,128,field) for field in fields}
summary_fieldnames=["Dataset","Target PSNR","Method","CPUs","Samples","Sampled Element Num","Compression Throughput MB/s","Decompression Throughput MB/s","Compression Throughput Per CPU MB/s","Decompression Throughput Per CPU MB/s","Compression Throughput Normalization Factor","Decompression Throughput Normalization Factor","Normalized Compression Throughput MB/s","Normalized Decompression Throughput MB/s","Compression Efficiency","Compression Efficiency Std Dev","Decompression Efficiency","Decompression Efficiency Std Dev","Aggregate BR"]
quality_fieldnames=["Dataset","Target PSNR","Method","Field","CPUs","REL EB from full RD","ABS EB","Aggregate BR","Aggregate MSE","PSNR using full range","Max Error"]

def read_rows(name):
    with open(os.path.join(result_dir,name),newline="") as f:
        return list(csv.DictReader(f))

def write_rows(name,fieldnames,rows):
    with open(os.path.join(result_dir,name),"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def row_identity(row):
    return int(row["Global Partition Begin"]),int(row["Global Partition End"]),row["Field"]

def matched_time_ratio(source_rows,reference_rows,column):
    source={row_identity(row):row for row in source_rows}
    matched=[(source[row_identity(row)],row) for row in reference_rows]
    if len(matched)!=len(reference_rows):
        raise RuntimeError("节点校正没有找到全部匹配case")
    return sum(float(left[column]) for left,_ in matched)/sum(float(right[column]) for _,right in matched)

def effective_time_factor(rows,column,node_factors):
    measured_time=sum(float(row[column]) for row in rows)
    normalized_time=sum(float(row[column])/node_factors["a245" if row_identity(row) in one_resume_identities else "a120"] for row in rows)
    return measured_time/normalized_time

def aggregate_throughput(rows,column):
    return sum(int(row["Raw Bytes"]) for row in rows)/sum(float(row[column]) for row in rows)/1024/1024

def validate(rows):
    for cpus,count in expected_rows.items():
        group=[row for row in rows if int(row["CPUs"])==cpus]
        if len(group)!=count:
            raise RuntimeError(f"{cpus} CPU应有{count}行，实际为{len(group)}行")
        identities={(row["Field"],row["Trial"],row["Element Begin"],row["Element End"]) for row in group}
        if len(identities)!=count:
            raise RuntimeError(f"{cpus} CPU存在重复case")
        covered=sum(int(row["Element Num"]) for row in group)//len(fields)
        expected=element_num*expected_repeats[cpus]
        if covered!=expected:
            raise RuntimeError(f"{cpus} CPU覆盖{covered}个element，预期{expected}")
        if any(float(row["Max Error"])>float(row["ABS EB"])*1.0001 for row in group):
            raise RuntimeError(f"{cpus} CPU存在误差越界")

def make_summary(rows,factors):
    baseline_rows=[row for row in rows if int(row["CPUs"])==1]
    baseline_compression=aggregate_throughput(baseline_rows,"Compression Time")*factors[1][0]
    baseline_decompression=aggregate_throughput(baseline_rows,"Decompression Time")*factors[1][1]
    summary=[]
    quality=[]
    for cpus in sorted(expected_rows):
        group=[row for row in rows if int(row["CPUs"])==cpus]
        compression=aggregate_throughput(group,"Compression Time")
        decompression=aggregate_throughput(group,"Decompression Time")
        compression_factor,decompression_factor=factors[cpus]
        normalized_compression=compression*compression_factor
        normalized_decompression=decompression*decompression_factor
        trial_compression=[]
        trial_decompression=[]
        for trial in sorted({int(row["Trial"]) for row in group}):
            trial_group=[row for row in group if int(row["Trial"])==trial]
            trial_compression.append(aggregate_throughput(trial_group,"Compression Time")*compression_factor/cpus/baseline_compression)
            trial_decompression.append(aggregate_throughput(trial_group,"Decompression Time")*decompression_factor/cpus/baseline_decompression)
        value_num=sum(int(row["Value Num"]) for row in group)
        summary.append({"Dataset":"jicf","Target PSNR":150.0,"Method":"FDWZ","CPUs":cpus,"Samples":len(group),"Sampled Element Num":sum(int(row["Element Num"]) for row in group),"Compression Throughput MB/s":compression,"Decompression Throughput MB/s":decompression,"Compression Throughput Per CPU MB/s":compression/cpus,"Decompression Throughput Per CPU MB/s":decompression/cpus,"Compression Throughput Normalization Factor":compression_factor,"Decompression Throughput Normalization Factor":decompression_factor,"Normalized Compression Throughput MB/s":normalized_compression,"Normalized Decompression Throughput MB/s":normalized_decompression,"Compression Efficiency":normalized_compression/cpus/baseline_compression,"Compression Efficiency Std Dev":statistics.stdev(trial_compression) if len(trial_compression)>1 else 0,"Decompression Efficiency":normalized_decompression/cpus/baseline_decompression,"Decompression Efficiency Std Dev":statistics.stdev(trial_decompression) if len(trial_decompression)>1 else 0,"Aggregate BR":sum(float(row["BR"])*int(row["Value Num"]) for row in group)/value_num})
        for field in fields:
            field_group=[row for row in group if row["Field"]==field]
            field_values=sum(int(row["Value Num"]) for row in field_group)
            mse=sum(float(row["MSE"])*int(row["Value Num"]) for row in field_group)/field_values
            quality.append({"Dataset":"jicf","Target PSNR":150.0,"Method":"FDWZ","Field":field,"CPUs":cpus,"REL EB from full RD":field_group[0]["REL EB from full RD"],"ABS EB":field_group[0]["ABS EB"],"Aggregate BR":sum(float(row["BR"])*int(row["Value Num"]) for row in field_group)/field_values,"Aggregate MSE":mse,"PSNR using full range":20*math.log10(ranges[field])-10*math.log10(mse),"Max Error":max(float(row["Max Error"]) for row in field_group)})
    return summary,quality

def normalized_existing_summary(rows):
    result=[]
    for row in rows:
        if row["Method"]=="FDWZ":
            continue
        row=dict(row)
        row["Compression Throughput Normalization Factor"]=1
        row["Decompression Throughput Normalization Factor"]=1
        row["Normalized Compression Throughput MB/s"]=row["Compression Throughput MB/s"]
        row["Normalized Decompression Throughput MB/s"]=row["Decompression Throughput MB/s"]
        result.append(row)
    return result

def main():
    one=read_rows("weak_jicf_current_fdwz_1cpu128_full_raw.csv")
    low=read_rows("weak_jicf_current_fdwz_full_2to16_raw.csv")
    reference_name="weak_jicf_current_fdwz_sampled_reference_raw.csv"
    reference_path=os.path.join(result_dir,reference_name)
    if os.path.isfile(reference_path):
        reference=read_rows(reference_name)
    else:
        reference=read_rows("weak_jicf_current_fdwz_raw.csv")
        write_rows(reference_name,list(reference[0]),reference)
    if len(reference)!=160:
        raise RuntimeError(f"原始同节点校正结果应有160行，实际为{len(reference)}行")
    high=[row for row in reference if int(row["CPUs"])>=32]
    merged=one+low+high
    validate(merged)
    reference_one=[row for row in reference if int(row["CPUs"])==1]
    reference_sixteen=[row for row in reference if int(row["CPUs"])==16]
    low_sixteen=[row for row in low if int(row["CPUs"])==16]
    reference_one_a120=[row for row in reference_one if row_identity(row) not in one_resume_identities]
    reference_one_a245=[row for row in reference_one if row_identity(row) in one_resume_identities]
    one_node_factors={}
    for column in ["Compression Time","Decompression Time"]:
        one_node_factors[column]={"a120":matched_time_ratio(one,reference_one_a120,column),"a245":matched_time_ratio(one,reference_one_a245,column)}
    one_factors=tuple(effective_time_factor(one,column,one_node_factors[column]) for column in ["Compression Time","Decompression Time"])
    low_factors=(matched_time_ratio(low_sixteen,reference_sixteen,"Compression Time"),matched_time_ratio(low_sixteen,reference_sixteen,"Decompression Time"))
    factors={1:one_factors,2:low_factors,4:low_factors,8:low_factors,16:low_factors,32:(1,1),64:(1,1),128:(1,1)}
    summary,quality=make_summary(merged,factors)
    raw_fieldnames=list(merged[0])
    write_rows("weak_jicf_current_fdwz_raw.csv",raw_fieldnames,merged)
    write_rows("weak_jicf_current_fdwz_summary.csv",summary_fieldnames,summary)
    write_rows("weak_jicf_current_fdwz_quality.csv",quality_fieldnames,quality)
    calibration=[{"Source Node":"a120","Reference Node":"a217","Matched CPUs":1,"Matched Cases":len(reference_one_a120),"Compression Throughput Normalization Factor":one_node_factors["Compression Time"]["a120"],"Decompression Throughput Normalization Factor":one_node_factors["Decompression Time"]["a120"]},{"Source Node":"a245","Reference Node":"a217","Matched CPUs":1,"Matched Cases":len(reference_one_a245),"Compression Throughput Normalization Factor":one_node_factors["Compression Time"]["a245"],"Decompression Throughput Normalization Factor":one_node_factors["Decompression Time"]["a245"]},{"Source Node":"a124","Reference Node":"a217","Matched CPUs":16,"Matched Cases":len(reference_sixteen),"Compression Throughput Normalization Factor":low_factors[0],"Decompression Throughput Normalization Factor":low_factors[1]}]
    write_rows("weak_jicf_current_fdwz_node_calibration.csv",list(calibration[0]),calibration)
    formal_raw=[row for row in read_rows("weak_jicf_raw.csv") if row["Method"]!="FDWZ"]+merged
    method_order={method:index for index,method in enumerate(methods)}
    field_order={field:index for index,field in enumerate(fields)}
    formal_raw.sort(key=lambda row:(method_order[row["Method"]],int(row["CPUs"]),int(row["Trial"]),field_order.get(row["Field"],0)))
    write_rows("weak_jicf_raw.csv",raw_fieldnames,formal_raw)
    formal_summary=normalized_existing_summary(read_rows("weak_jicf_summary.csv"))+summary
    formal_summary.sort(key=lambda row:(method_order[row["Method"]],int(row["CPUs"])))
    write_rows("weak_jicf_summary.csv",summary_fieldnames,formal_summary)
    formal_quality=[row for row in read_rows("weak_jicf_quality.csv") if row["Method"]!="FDWZ"]+quality
    formal_quality.sort(key=lambda row:(method_order[row["Method"]],int(row["CPUs"]),field_order.get(row["Field"],0)))
    write_rows("weak_jicf_quality.csv",quality_fieldnames,formal_quality)
    print(f"merged {len(merged)} FDWZ rows")
    print(f"normalization factors: one_effective={one_factors}, one_nodes={one_node_factors}, low={low_factors}")

if __name__=="__main__":
    main()
