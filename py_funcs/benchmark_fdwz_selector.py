import argparse,csv,glob,os,re,subprocess,sys
from pathlib import Path

project_root=Path(__file__).resolve().parents[1]
ansi_re=re.compile(r"\x1b\[[0-9;]*m")

dataset_configs={
    "pb0":("nekrs",Path("/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001"),None,["UX","UY","UZ","P","T"]),
    "mixlay0":("nekrs",project_root/"third_party/pySEMTools/examples/data/mixlay0.f00001",None,["UX","UY","P","T","S00","S01"]),
    "test0":("nekrs",project_root/"third_party/fld_pyutils/demos/data/test0.f00001",None,["UZ","T","S00","S01"]),
    "tc_channel0":("nekrs",project_root/"third_party/pySEMTools/examples/data/tc_channel0.f00001",None,["UX","UY","UZ"]),
    "mfem_ex0":("mfem",project_root/"build/mfem_build/examples/ex0examples/sol.gf",project_root/"build/mfem_build/examples/ex0examples/mesh.mesh",None),
    "mfem_ex1":("mfem",project_root/"build/mfem_build/examples/ex1examples/sol.gf",project_root/"build/mfem_build/examples/ex1examples/refined.mesh",None),
    "nektar_cylinder":("nektar++",project_root/"third_party/flow-stability-cylinder/tutorial/base/Cylinder-Base.vtu",None,["p","u","v"]),
}

def clean_line(line):
    return ansi_re.sub("",line).strip()

def run_command(command):
    result=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
    print(result.stdout,end="",flush=True)
    return result.stdout

def parse_compression(output):
    results={}
    current_field=None
    for raw_line in output.splitlines():
        line=clean_line(raw_line)
        if line.startswith("Field= ") and "," not in line:
            current_field=line.split("=",1)[1].strip()
            results.setdefault(current_field,{})
        elif current_field is not None and line.startswith("Sampled configurations:"):
            results[current_field]["Configuration"]=line.split(":",1)[1].strip()
        elif current_field is not None and line.startswith("FDWZ QB width= "):
            results[current_field]["QB Width"]=int(line.split("=",1)[1].split()[0])
        elif current_field is not None and line.startswith("Compression time= "):
            results[current_field]["Compression Time"]=float(line.split("=",1)[1].split()[0])
        elif line.startswith("Field= ") and "compression bandwidth per cpu=" in line:
            field=line.split(",",1)[0].split("=",1)[1].strip()
            results.setdefault(field,{})["Compression Throughput"]=float(line.split("compression bandwidth per cpu=",1)[1].split()[0])
    return results

def parse_decompression(output):
    results={}
    for raw_line in output.splitlines():
        line=clean_line(raw_line)
        if line.startswith("Field= ") and "decompression bandwidth per cpu=" in line:
            field=line.split(",",1)[0].split("=",1)[1].strip()
            results.setdefault(field,{})["Decompression Throughput"]=float(line.split("decompression bandwidth per cpu=",1)[1].split()[0])
        elif line.startswith("Field= ") and ", abs_eb=" in line:
            segments=[segment.strip() for segment in line.split(",")]
            field=segments[0].split("=",1)[1].strip()
            values=results.setdefault(field,{})
            values["BR"]=float(segments[3].split("=",1)[1])
            values["CR"]=float(segments[4].split("=",1)[1])
            values["Max Error Percent"]=float(segments[5].rsplit("(",1)[1].removesuffix("%)"))
            values["PSNR"]=float(segments[7].split("=",1)[1].split()[0])
    return results

def cleanup_outputs(data_path,method):
    for path in glob.glob(f"{glob.escape(str(data_path))}_*.{glob.escape(method)}"):
        Path(path).unlink(missing_ok=True)
    for suffix in [f".{method}.bin","_header.mrcz","_element_map.mrcz"]:
        Path(f"{data_path}{suffix}").unlink(missing_ok=True)

def build_common_command(binary,format_name,data_path,mesh_path,fields,eb,method,max_cpus):
    command=[str(project_root/"install/bin"/binary),"-format",format_name,"-i",str(data_path)]
    if mesh_path is not None:
        command.extend(["-mesh",str(mesh_path)])
    command.extend(["-M","REL",str(eb),"-method",method,"-max-cpus",str(max_cpus),"-partition-num",str(max_cpus)])
    if fields is not None:
        command.extend(["-fields",f"[{','.join(fields)}]"])
    return command

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--dataset",choices=dataset_configs,required=True)
    parser.add_argument("--ebs",default="1e-4,1e-5")
    parser.add_argument("--fields")
    parser.add_argument("--max-cpus",type=int,default=1)
    parser.add_argument("--repeats",type=int,default=1)
    parser.add_argument("--sample-interval",type=int,default=128)
    parser.add_argument("--measure-phase-time",action="store_true")
    parser.add_argument("--strategies",default="linear-proxy,cubic-exhaustive")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.max_cpus<1 or args.repeats<1 or args.sample_interval<1:
        parser.error("--max-cpus、--repeats和--sample-interval必须为正整数")
    format_name,data_path,mesh_path,fields=dataset_configs[args.dataset]
    if args.fields is not None:
        fields=[field for field in args.fields.split(",") if field]
    ebs=[float(value) for value in args.ebs.split(",") if value]
    strategies=[value for value in args.strategies.split(",") if value]
    fixed_plans={f"{order}-cubic" for order in ["XYZ","XZY","YXZ","YZX","ZXY","ZYX"]}
    allowed={"linear-only","linear-proxy","cubic-exhaustive","SZ3-4D"}|fixed_plans
    if any(strategy not in allowed for strategy in strategies):
        parser.error("--strategies包含不支持的策略")
    rows=[]
    for eb in ebs:
        for repeat in range(args.repeats):
            ordered_strategies=strategies[repeat%len(strategies):]+strategies[:repeat%len(strategies)]
            for strategy in ordered_strategies:
                method="SZ3-4D" if strategy=="SZ3-4D" else "FDWZ-ablation" if strategy in fixed_plans else "FDWZ"
                print(f"BENCHMARK dataset={args.dataset} eb={eb:g} repeat={repeat+1} strategy={strategy}",flush=True)
                compression_command=build_common_command("mrcz_compress",format_name,data_path,mesh_path,fields,eb,method,args.max_cpus)
                if args.measure_phase_time:
                    compression_command.extend(["-measure-phase-time","true"])
                if method=="FDWZ":
                    compression_command.extend(["-fdwz-sample-interval",str(args.sample_interval),"-fdwz-selector",strategy])
                elif method=="FDWZ-ablation":
                    compression_command.extend(["-fdwz-plan",strategy])
                decompression_command=build_common_command("mrcz_decompress",format_name,data_path,mesh_path,fields,eb,method,args.max_cpus)
                try:
                    compression=parse_compression(run_command(compression_command))
                    decompression=parse_decompression(run_command(decompression_command))
                    for field in sorted(compression):
                        values=compression[field]|decompression.get(field,{})
                        rows.append({"Dataset":args.dataset,"REL EB":eb,"Repeat":repeat+1,"Strategy":strategy,"Field":field,**values})
                finally:
                    cleanup_outputs(data_path,method)
    columns=["Dataset","REL EB","Repeat","Strategy","Field","Configuration","QB Width","Compression Time","Compression Throughput","Decompression Throughput","BR","CR","PSNR","Max Error Percent"]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(args.output.resolve())

if __name__=="__main__":
    main()
