import csv,os,sys

project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,project_root)
from t5_start_weak_scaling import datasets
from t6_start_strong_scaling import fieldnames,normalize_strong_row,parallel_methods,write_rows,write_summaries

strong_dir=os.path.join(project_root,"results","strong_scaling")
weak_dir=os.path.join(project_root,"results","weak_scaling")
source_name="strong_pb0_8to32_raw.csv"
formal_name="strong_pb0_raw.csv"

def read_rows(path):
    with open(path,newline="") as file:
        return list(csv.DictReader(file))

def main():
    source_path=os.path.join(strong_dir,source_name)
    formal_path=os.path.join(strong_dir,formal_name)
    if not os.path.isfile(source_path):
        source_rows=read_rows(formal_path)
        source_rows=[row for row in source_rows if int(row["CPUs"])<64]
        if len(source_rows)!=90:
            raise RuntimeError(f"8--32 CPU strong结果应有90行，实际为{len(source_rows)}行")
        write_rows(source_path,list(source_rows[0]),source_rows)
    strong_rows=[normalize_strong_row(row) for row in read_rows(source_path)]
    for row in strong_rows:
        row["Source Job"]="20454976"
        row["Source Node"]="a095"
        row["Source Result"]=source_name
    weak_rows=read_rows(os.path.join(weak_dir,"weak_pb0_raw.csv"))
    weak64=[row for row in weak_rows if row["Method"] in parallel_methods and int(row["CPUs"])==64]
    if len(weak64)!=120:
        raise RuntimeError(f"64 CPU weak复用结果应有120行，实际为{len(weak64)}行")
    imported=[]
    config=datasets["pb0"]
    for row in weak64:
        if int(row["Element Begin"])!=0 or int(row["Element End"])!=config["element_num"] or int(row["Partition Num"])!=64:
            raise RuntimeError("64 CPU weak结果不是完整数据和固定64 partition")
        result={column:row[column] for column in fieldnames if column in row}
        result["Source Job"]="20424829" if row["Method"]=="FDWZ" else "20348291"
        result["Source Node"]="a173" if row["Method"]=="FDWZ" else "a214"
        result["Source Result"]="weak_pb0_raw.csv"
        imported.append(normalize_strong_row(result))
    merged=strong_rows+imported
    method_order={method:index for index,method in enumerate(parallel_methods)}
    field_order={field:index for index,field in enumerate(config["fields"])}
    merged.sort(key=lambda row:(method_order[row["Method"]],row["CPUs"],row["Trial"],field_order[row["Field"]]))
    for method in parallel_methods:
        for field in config["fields"]:
            group=[row for row in merged if row["Method"]==method and row["Field"]==field]
            if len(group)!=7:
                raise RuntimeError(f"{method}/{field}应有3个新CPU点和4次64 CPU复用结果")
            if len({row["ABS EB"] for row in group})!=1 or len({row["Compute Precision"] for row in group})!=1 or len({row["REL EB from full RD"] for row in group})!=1:
                raise RuntimeError(f"{method}/{field}在不同CPU数下的error bound或计算精度不一致")
            if max(row["Max Error"]/row["ABS EB"] for row in group)>1.0001:
                raise RuntimeError(f"{method}/{field}存在误差越界")
            if max(row["MSE"] for row in group)!=min(row["MSE"] for row in group):
                raise RuntimeError(f"{method}/{field}在不同CPU数下MSE不一致")
            if max(row["BR"] for row in group)-min(row["BR"] for row in group)>1.1e-6:
                raise RuntimeError(f"{method}/{field}在不同CPU数下BR不一致")
    write_rows(formal_path,fieldnames,merged)
    write_summaries(merged,config,"pb0","")
    print(f"merged {len(strong_rows)} strong rows and {len(imported)} reused 64-CPU rows")

if __name__=="__main__":
    main()
