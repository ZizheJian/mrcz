import csv,os,re,statistics,subprocess

project_root=os.path.abspath(os.path.join(os.path.dirname(__file__),"../.."))
compress_path=os.path.join(project_root,"install/bin/mrcz_compress")
decompress_path=os.path.join(project_root,"install/bin/mrcz_decompress")
output_dir=os.path.dirname(__file__)
raw_path=os.path.join(output_dir,"raw_results.csv")
summary_path=os.path.join(output_dir,"summary.csv")
log_path=os.path.join(output_dir,"raw.log")

methods=["MRCZ-3D-GLL","MRCZ-3D-GLL-coeff-delta-constant","MRCZ-3D-GLL-coeff-delta-linear","MRCZ-3D-GLL-coeff-delta-all"]
cases=[
    ("mixlay0",os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),["UX","UY","UZ","P","T"]),
    ("test0",os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),["T","S00","S01"]),
    ("tc_channel0",os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),["UX","UY","UZ"])]
eb_list=[1e-2,1e-3,1e-4,1e-5]
columns=["Dataset","Field","REL EB","Method","BR","PSNR","Max Error","ABS EB","Compression Time","Decompression Time"]

def run(command):
    result=subprocess.run(command,cwd=project_root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    with open(log_path,"a") as f:
        f.write(" ".join(command)+"\n"+result.stdout+"\n")
    if result.returncode!=0:
        raise subprocess.CalledProcessError(result.returncode,command,result.stdout)
    return result.stdout

def find_float(output,pattern):
    match=re.search(pattern,output)
    if match is None:
        raise RuntimeError("无法解析输出: "+pattern)
    return float(match.group(1))

def cleanup(data_path,method,field):
    for path in [f"{data_path}_{field}.{method}",f"{data_path}_header.mrcz",f"{data_path}_element_map.mrcz",f"{data_path}.{method}.bin"]:
        if os.path.isfile(path):
            os.remove(path)

def write_rows(rows):
    with open(raw_path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

def summarize(rows):
    baseline={(row["Dataset"],row["Field"],row["REL EB"]):row for row in rows if row["Method"]=="MRCZ-3D-GLL"}
    summary=[]
    for method in methods[1:]:
        selected=[row for row in rows if row["Method"]==method]
        ratios=[]
        psnr_differences=[]
        max_error_ratios=[]
        for row in selected:
            base=baseline[(row["Dataset"],row["Field"],row["REL EB"])]
            ratios.append(row["BR"]/base["BR"])
            psnr_differences.append(abs(row["PSNR"]-base["PSNR"]))
            max_error_ratios.append(row["Max Error"]/row["ABS EB"])
        summary.append({"Method":method,"Cases":len(selected),"BR Geometric Mean Change (%)":(statistics.geometric_mean(ratios)-1)*100,"BR Wins":sum(ratio<1 for ratio in ratios),"BR Losses":sum(ratio>1 for ratio in ratios),"Best BR Change (%)":(min(ratios)-1)*100,"Worst BR Change (%)":(max(ratios)-1)*100,"Max PSNR Difference":max(psnr_differences),"Max Error/EB":max(max_error_ratios)})
    with open(summary_path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

rows=[]
if os.path.exists(log_path):
    os.remove(log_path)
for dataset,data_path,fields in cases:
    for field in fields:
        for eb in eb_list:
            for method in methods:
                common=["-format","nekrs","-i",data_path,"-M","REL",str(eb),"-method",method,"-max-cpus","1","-partition-num","1","-fields",f"[{field}]"]
                compression_output=run([compress_path]+common)
                decompression_output=run([decompress_path]+common)
                rows.append({"Dataset":dataset,"Field":field,"REL EB":eb,"Method":method,"BR":find_float(decompression_output,r"BR= ([0-9eE+.-]+)"),"PSNR":find_float(decompression_output,r"PSNR= ([0-9eE+.-]+)"),"Max Error":find_float(decompression_output,r"max_err= ([0-9eE+.-]+)"),"ABS EB":find_float(decompression_output,r"abs_eb= ([0-9eE+.-]+)"),"Compression Time":find_float(compression_output,r"Compression time= ([0-9eE+.-]+)"),"Decompression Time":find_float(decompression_output,r"Decompression time= ([0-9eE+.-]+)")})
                write_rows(rows)
                cleanup(data_path,method,field)
                print(f"完成: {dataset} {field} REL={eb} {method}")
summarize(rows)
