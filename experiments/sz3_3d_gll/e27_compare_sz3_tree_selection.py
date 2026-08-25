import csv,os,re,shutil,tempfile,time
import numpy as np
from py_funcs.read_nekrs_header import read_nekrs_header
from third_party.run_command.run_command import run_command

# XXX 在更多field和error bound上比较SZ3-3D-GLL的Huffman tree选择方法。
# XXX 每个case先用分散sample选定维度顺序和linear/cubic，之后固定这两项，只计算tree数量导致的holdout BR regret。
# XXX 比较直接使用sample BR、用半sample和完sample的两点线性外推、固定2tree、固定最大tree和两种符号数阈值规则。
# XXX 两点外推将bitstream bits拟合为bits=per_element_bits*element_num+fixed_metadata_bits，再外推到完整partition。
# XXX 所有sample、half sample和holdout使用同一partition数值范围确定的ABS error bound，且所有临时输入和输出都位于/tmp。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
orders=["XYZ","XZY","YXZ","YZX","ZXY","ZYX"]
datasets=[
    {"name":"test0","path":os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),"fields":["T"],"rel_ebs":[1e-3,1e-4,1e-5],"max_trees":4,"window":(0,512,128)},
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"fields":["T"],"rel_ebs":[1e-3,1e-4,1e-5],"max_trees":5,"window":(0,1600,400)},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"fields":["UZ"],"rel_ebs":[1e-3,1e-4,1e-5],"max_trees":5,"window":(0,512,128)},
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","fields":["UZ","P","T"],"rel_ebs":[1e-3,1e-4,1e-5],"max_trees":5,"window":(1335789,1353197,1024)},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","fields":["UZ","P","S00"],"rel_ebs":[1e-4,1e-6,1e-8],"max_trees":5,"window":(578353,595761,1024)},
]
old_symbol_threshold=106496
natural_symbol_threshold=1000000
detail_path=os.path.join(experiment_dir,"e27_sz3_tree_selection_details.csv")
result_path=os.path.join(experiment_dir,"e27_sz3_tree_selection_results.csv")
summary_path=os.path.join(experiment_dir,"e27_sz3_tree_selection_summary.csv")

def method_name(order,predictor,tree_num):
    return f"SZ3-3D-GLL-sampled-{order}-{predictor}-{tree_num}"

def parse_number(output,pattern,name):
    match=re.search(pattern,output)
    if match is None:
        raise RuntimeError(f"无法从输出中解析{name}")
    return float(match.group(1))

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup_outputs(data_path,method):
    remove_if_exists(f"{data_path}_T.{method}")
    remove_if_exists(f"{data_path}.{method}.bin")

def open_field(path,field):
    with open(path,"rb") as source_file:
        header=read_nekrs_header(source_file)
        source_file.seek(0)
        original_header=source_file.read(136)
    dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
    offset=136+header.file_element_num*4
    field_value_num=header.file_element_num*header.element_size
    for variable in header.variables:
        if variable=="X":
            offset+=3*field_value_num*header.word_size
        elif variable=="U":
            if field in ["UX","UY","UZ"]:
                component=["UX","UY","UZ"].index(field)
                values=np.memmap(path,dtype=dtype,mode="r",offset=offset,shape=(header.file_element_num,3,header.element_size))
                return header,original_header,values[:,component,:]
            offset+=3*field_value_num*header.word_size
        elif variable in ["P","T"]:
            if field==variable:
                values=np.memmap(path,dtype=dtype,mode="r",offset=offset,shape=(header.file_element_num,header.element_size))
                return header,original_header,values
            offset+=field_value_num*header.word_size
        elif variable=="S":
            for scalar_id in range(header.passive_scalar_num):
                scalar_name=f"S{scalar_id:02d}"
                if field==scalar_name:
                    values=np.memmap(path,dtype=dtype,mode="r",offset=offset,shape=(header.file_element_num,header.element_size))
                    return header,original_header,values
                offset+=field_value_num*header.word_size
    raise ValueError(f"{path}不包含field {field}")

def write_compact_nekrs(path,header,original_header,values,indices):
    header_terms=original_header[:132].split()
    header_terms[5]=str(len(indices)).encode()
    header_terms[6]=str(len(indices)).encode()
    header_terms[11]=b"T"
    new_header=b" ".join(header_terms)
    if len(new_header)>132:
        raise RuntimeError("重写后的NekRS header超过132 bytes")
    map_dtype=np.dtype(np.int32).newbyteorder(header.endian_mode)
    with open(path,"wb") as output_file:
        output_file.write(new_header.ljust(132,b" "))
        output_file.write(original_header[132:136])
        np.arange(1,len(indices)+1,dtype=map_dtype).tofile(output_file)
        np.asarray(values[indices]).tofile(output_file)

def extrapolate_bits(half_bits,half_num,sample_bits,sample_num,partition_num):
    per_element_bits=(sample_bits-half_bits)/(sample_num-half_num)
    metadata_bits=sample_bits-per_element_bits*sample_num
    if per_element_bits<0 or metadata_bits<0:
        return sample_bits*partition_num/sample_num
    return per_element_bits*partition_num+metadata_bits

detail_fields=["Dataset","Field","REL EB","Order","Predictor","Tree","Half-sample BR","Sample BR","Holdout BR","Estimated Partition Bits"]
result_fields=["Dataset","Field","REL EB","Order","Predictor","Best Tree","Policy","Selected Tree","Best Holdout BR","Selected Holdout BR","Regret (%)"]
policies=["sample","two-point extrapolation","fixed 2tree","fixed maximum","old symbol threshold","1M symbol threshold"]
detail_rows=[]
result_rows=[]
temp_dir=tempfile.mkdtemp(prefix="mrcz-e27-")
experiment_start=time.perf_counter()
try:
    for dataset_id,dataset in enumerate(datasets):
        partition_begin,partition_end,sample_num=dataset["window"]
        partition_num=partition_end-partition_begin
        relative_indices=np.arange(partition_num,dtype=np.int64)
        rng=np.random.default_rng(2700+dataset_id)
        random_sample_indices=rng.choice(relative_indices,size=sample_num,replace=False)
        half_num=sample_num//2
        half_indices=np.sort(random_sample_indices[:half_num]+partition_begin)
        sample_indices=np.sort(random_sample_indices+partition_begin)
        holdout_indices=np.setdiff1d(relative_indices,random_sample_indices,assume_unique=False)+partition_begin
        for field in dataset["fields"]:
            header,original_header,field_values=open_field(dataset["path"],field)
            role_paths={role:os.path.join(temp_dir,f"{dataset['name']}-{field}-{role}.f00001") for role in ["half","sample","holdout"]}
            write_compact_nekrs(role_paths["half"],header,original_header,field_values,half_indices)
            write_compact_nekrs(role_paths["sample"],header,original_header,field_values,sample_indices)
            write_compact_nekrs(role_paths["holdout"],header,original_header,field_values,holdout_indices)
            partition_range=float(np.max(field_values[partition_begin:partition_end])-np.min(field_values[partition_begin:partition_end]))
            for rel_eb in dataset["rel_ebs"]:
                abs_eb=partition_range*rel_eb
                cache={}
                def measure(role,order,predictor,tree_num):
                    key=(role,order,predictor,tree_num)
                    if key in cache:
                        return cache[key]
                    data_path=role_paths[role]
                    method=method_name(order,predictor,tree_num)
                    cleanup_outputs(data_path,method)
                    args=f"-format nekrs -i {data_path} -M ABS {abs_eb:.17g} -method {method} -max-cpus 1 -fields [T]"
                    compression_output=run_command(f"./install/bin/mrcz_compress {args}",silent=True)
                    decompression_output=run_command(f"./install/bin/mrcz_decompress {args}",silent=True)
                    result={"BR":parse_number(decompression_output,r"BR= ([0-9.eE+-]+)","BR"),"PSNR":parse_number(decompression_output,r"PSNR= ([0-9.eE+-]+)","PSNR"),"Max Error":parse_number(decompression_output,r"max_err= ([0-9.eE+-]+)","max error")}
                    cache[key]=result
                    cleanup_outputs(data_path,method)
                    return result

                order=min(orders,key=lambda candidate:measure("sample",candidate,"linear",2)["BR"])
                predictor=min(["linear","cubic"],key=lambda candidate:measure("sample",order,candidate,2)["BR"])
                tree_results={}
                for tree_num in range(1,dataset["max_trees"]+1):
                    half_result=measure("half",order,predictor,tree_num)
                    sample_result=measure("sample",order,predictor,tree_num)
                    holdout_result=measure("holdout",order,predictor,tree_num)
                    half_bits=half_result["BR"]*half_num*header.element_size
                    sample_bits=sample_result["BR"]*sample_num*header.element_size
                    estimated_bits=extrapolate_bits(half_bits,half_num,sample_bits,sample_num,partition_num)
                    tree_results[tree_num]={"half":half_result,"sample":sample_result,"holdout":holdout_result,"estimated_bits":estimated_bits}
                    detail_rows.append({"Dataset":dataset["name"],"Field":field,"REL EB":rel_eb,"Order":order,"Predictor":predictor,"Tree":tree_num,"Half-sample BR":half_result["BR"],"Sample BR":sample_result["BR"],"Holdout BR":holdout_result["BR"],"Estimated Partition Bits":estimated_bits})
                best_tree=min(tree_results,key=lambda tree_num:tree_results[tree_num]["holdout"]["BR"])
                sample_tree=min(tree_results,key=lambda tree_num:tree_results[tree_num]["sample"]["BR"])
                extrapolated_tree=min(tree_results,key=lambda tree_num:tree_results[tree_num]["estimated_bits"])
                symbol_num=partition_num*header.element_size
                selected_trees={
                    "sample":sample_tree,"two-point extrapolation":extrapolated_tree,"fixed 2tree":2,"fixed maximum":dataset["max_trees"],
                    "old symbol threshold":dataset["max_trees"] if symbol_num>=old_symbol_threshold else 2,
                    "1M symbol threshold":dataset["max_trees"] if symbol_num>=natural_symbol_threshold else 2,
                }
                best_br=tree_results[best_tree]["holdout"]["BR"]
                for policy in policies:
                    selected_tree=selected_trees[policy]
                    selected_br=tree_results[selected_tree]["holdout"]["BR"]
                    result_rows.append({"Dataset":dataset["name"],"Field":field,"REL EB":rel_eb,"Order":order,"Predictor":predictor,"Best Tree":best_tree,"Policy":policy,"Selected Tree":selected_tree,"Best Holdout BR":best_br,"Selected Holdout BR":selected_br,"Regret (%)":(selected_br/best_br-1)*100})
                print(f"{dataset['name']} {field} REL={rel_eb:g}: {order}/{predictor},best={best_tree},sample={sample_tree},extrapolated={extrapolated_tree}",flush=True)
            for role_path in role_paths.values():
                remove_if_exists(f"{role_path}_header.mrcz")
                remove_if_exists(f"{role_path}_element_map.mrcz")
            del field_values
finally:
    shutil.rmtree(temp_dir)

with open(detail_path,"w",newline="") as output_file:
    writer=csv.DictWriter(output_file,fieldnames=detail_fields)
    writer.writeheader()
    writer.writerows(detail_rows)
with open(result_path,"w",newline="") as output_file:
    writer=csv.DictWriter(output_file,fieldnames=result_fields)
    writer.writeheader()
    writer.writerows(result_rows)

summary_fields=["Policy","Case Count","Exact Count","Mean Regret (%)","Max Regret (%)","Aggregate BR Regret (%)"]
summary_rows=[]
for policy in policies:
    policy_rows=[row for row in result_rows if row["Policy"]==policy]
    summary_rows.append({
        "Policy":policy,"Case Count":len(policy_rows),"Exact Count":sum(row["Selected Tree"]==row["Best Tree"] for row in policy_rows),
        "Mean Regret (%)":sum(row["Regret (%)"] for row in policy_rows)/len(policy_rows),"Max Regret (%)":max(row["Regret (%)"] for row in policy_rows),
        "Aggregate BR Regret (%)":(sum(row["Selected Holdout BR"] for row in policy_rows)/sum(row["Best Holdout BR"] for row in policy_rows)-1)*100,
    })
with open(summary_path,"w",newline="") as output_file:
    writer=csv.DictWriter(output_file,fieldnames=summary_fields)
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"实验用时{time.perf_counter()-experiment_start:.3f}秒")
for row in summary_rows:
    print(f"{row['Policy']}: exact={row['Exact Count']}/{row['Case Count']},mean regret={row['Mean Regret (%)']:.4f}%,max={row['Max Regret (%)']:.4f}%,aggregate={row['Aggregate BR Regret (%)']:.4f}%")
