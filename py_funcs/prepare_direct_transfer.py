import argparse,csv,hashlib,os

PROJECT_ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTPUT_ROOT="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/direct_transfer"
DATASETS={
    "pb0":{"source":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/hyuan/pb0.f00001","fields":["UX","UY","UZ","P","T"]},
    "jicf0":{"source":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/viral-flamespeed/jicf0.f00000","fields":["UX","UY","UZ","P","S00"]},
}

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--dataset",choices=["all",*DATASETS],default="all")
    parser.add_argument("--output-root",default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--overwrite",action="store_true")
    parser.add_argument("--dry-run",action="store_true")
    return parser.parse_args()

def read_header(path):
    with open(path,"rb") as input_file:header=input_file.read(136)
    if len(header)!=136:raise RuntimeError(f"{path}的NekRS header不足136字节")
    terms=header[:132].decode("ascii").split()
    if len(terms)<12 or terms[0]!="#std":raise RuntimeError(f"{path}的NekRS header格式无法识别")
    word_size=int(terms[1])
    shape=tuple(int(value) for value in terms[2:5])
    element_num=int(terms[6])
    variables=terms[11]
    dimension=3 if shape[2]>1 else 2 if shape[1]>1 else 1
    field_bytes=element_num*shape[0]*shape[1]*shape[2]*word_size
    offset=136+element_num*4
    fields={}
    variable_index=0
    while variable_index<len(variables):
        variable=variables[variable_index]
        if variable=="X":offset+=field_bytes*dimension
        elif variable=="U":
            for component,name in enumerate(("UX","UY","UZ")[:dimension]):fields[name]=(offset,component)
            offset+=field_bytes*dimension
        elif variable in ("P","T"):
            fields[variable]=(offset,None)
            offset+=field_bytes
        elif variable=="S":
            scalar_num=int(variables[variable_index+1:variable_index+3])
            for scalar_index in range(scalar_num):
                fields[f"S{scalar_index:02d}"]=(offset,None)
                offset+=field_bytes
            variable_index+=2
        variable_index+=1
    return {"field_bytes":field_bytes,"fields":fields}

def hash_file(path):
    digest=hashlib.sha256()
    with open(path,"rb") as input_file:
        while payload:=input_file.read(64*1024*1024):digest.update(payload)
    return digest.hexdigest()

def export_dataset(dataset_name,config,header,output_root,overwrite):
    fields=config["fields"]
    missing=[field for field in fields if field not in header["fields"]]
    if missing:raise RuntimeError(f"{config['source']}缺少field: {','.join(missing)}")
    first_offset=header["fields"][fields[0]][0]
    expected_offsets=[first_offset+index*header["field_bytes"] for index in range(len(fields))]
    actual_offsets=[header["fields"][field][0]+(header["fields"][field][1] or 0)*header["field_bytes"] for field in fields]
    if actual_offsets!=expected_offsets:raise RuntimeError(f"{config['source']}中的目标field不是连续存储")
    total_bytes=header["field_bytes"]*len(fields)
    output_directory=os.path.join(output_root,dataset_name)
    os.makedirs(output_directory,exist_ok=True)
    output_path=os.path.join(output_directory,f"{dataset_name}.bin")
    if os.path.exists(output_path) and not overwrite:
        if os.path.getsize(output_path)!=total_bytes:raise RuntimeError(f"{output_path}已存在但大小不正确，请使用--overwrite")
        return output_path,total_bytes,hash_file(output_path)
    temporary_path=output_path+".part"
    if os.path.exists(temporary_path):os.remove(temporary_path)
    digest=hashlib.sha256()
    remaining=total_bytes
    with open(config["source"],"rb") as input_file,open(temporary_path,"wb") as output_file:
        input_file.seek(first_offset)
        while remaining:
            payload=input_file.read(min(64*1024*1024,remaining))
            if not payload:raise RuntimeError(f"读取{config['source']}时提前到达文件末尾")
            output_file.write(payload)
            digest.update(payload)
            remaining-=len(payload)
    if os.path.getsize(temporary_path)!=total_bytes:raise RuntimeError(f"{temporary_path}大小不正确")
    os.replace(temporary_path,output_path)
    return output_path,total_bytes,digest.hexdigest()

def write_manifest(manifest):
    manifest_path=os.path.join(PROJECT_ROOT,"results","end_to_end","direct_transfer_manifest.csv")
    with open(manifest_path,"w",newline="") as output_file:
        writer=csv.DictWriter(output_file,fieldnames=["Dataset","Fields","Path","Bytes","SHA256"],lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    print(manifest_path)

def main():
    args=parse_args()
    selected=list(DATASETS) if args.dataset=="all" else [args.dataset]
    manifest=[]
    for dataset_name in selected:
        config=DATASETS[dataset_name]
        header=read_header(config["source"])
        total_bytes=header["field_bytes"]*len(config["fields"])
        print(f"{dataset_name}: 1 file containing {len(config['fields'])} fields, {total_bytes} bytes")
        if args.dry_run:continue
        output_path,total_bytes,digest=export_dataset(dataset_name,config,header,args.output_root,args.overwrite)
        manifest.append({"Dataset":dataset_name,"Fields":";".join(config["fields"]),"Path":output_path,"Bytes":total_bytes,"SHA256":digest})
        print(f"{output_path}: sha256={digest}")
    if not args.dry_run:write_manifest(manifest)

if __name__=="__main__":main()
