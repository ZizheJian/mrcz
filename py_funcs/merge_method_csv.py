import argparse,csv,io,subprocess
from pathlib import Path

project_root=Path(__file__).resolve().parents[1]

def read_csv(text):
    reader=csv.DictReader(io.StringIO(text))
    return list(reader.fieldnames or []),list(reader)

def row_key(row,key_columns):
    return tuple(float(row[column]) for column in key_columns)

def read_base(destination,base_ref):
    if base_ref is None:
        return destination.read_text()
    relative_path=destination.resolve().relative_to(project_root)
    return subprocess.run(["git","show",f"{base_ref}:{relative_path}"],cwd=project_root,text=True,stdout=subprocess.PIPE,check=True).stdout

def merge_file(source,destination,method,base_ref,allow_partial):
    source_columns,source_rows=read_csv(source.read_text())
    base_columns,base_rows=read_csv(read_base(destination,base_ref))
    method_columns=[column for column in source_columns if column.startswith(f"{method}_") or column.startswith(f"{method} ")]
    if not method_columns:
        raise RuntimeError(f"{source}中没有{method}结果列")
    key_columns=["REL EB"]
    if "Max CPUs" in source_columns and "Max CPUs" in base_columns:
        key_columns.append("Max CPUs")
    if "Partition Num" in source_columns and "Partition Num" in base_columns:
        key_columns.append("Partition Num")
    source_by_key={row_key(row,key_columns):row for row in source_rows}
    updated=0
    for row in base_rows:
        source_row=source_by_key.get(row_key(row,key_columns))
        if source_row is None:
            continue
        for column in method_columns:
            row[column]=source_row[column]
            if column not in base_columns:
                base_columns.append(column)
        updated+=1
    if updated==0 or not allow_partial and updated!=len(source_rows):
        raise RuntimeError(f"{source}只有{updated}/{len(source_rows)}行匹配{destination}")
    with destination.open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=base_columns,lineterminator="\n")
        writer.writeheader()
        writer.writerows(base_rows)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("source_directory",type=Path)
    parser.add_argument("destination_directory",type=Path)
    parser.add_argument("--method",required=True)
    parser.add_argument("--pattern",default="*.csv")
    parser.add_argument("--base-ref")
    parser.add_argument("--allow-partial",action="store_true")
    args=parser.parse_args()
    sources=sorted(args.source_directory.glob(args.pattern))
    if not sources:
        parser.error("没有匹配的source CSV")
    for source in sources:
        destination=args.destination_directory/source.name
        if not destination.is_file():
            raise RuntimeError(f"缺少目标CSV: {destination}")
        merge_file(source,destination,args.method,args.base_ref,args.allow_partial)
        print(destination.resolve())

if __name__=="__main__":
    main()
