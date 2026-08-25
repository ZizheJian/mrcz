import csv,os,re,shutil,tempfile
import numpy as np
from py_funcs.read_nekrs_header import read_nekrs_header
from third_party.run_command.run_command import run_command

# XXX 检查少量element采样能否为同一局部partition选择SZ3-3D-GLL配置。
# XXX 严格按照维度顺序、预测函数、Huffman tree数量的顺序选择：先在linear+2tree下比较六种顺序，再比较linear/cubic，最后比较1到4或5tree。
# XXX 从每个连续partition窗口中用固定随机种子分散抽取sample，其余element作为holdout；两者分别独立选择。
# XXX 然后把sample选中的配置用于holdout，报告相对正式XYZ+linear+2tree的收益和相对holdout层次最优配置的regret。
# XXX 所有候选在同一REL error bound下按真实bitstream BR选择；PSNR和最大误差只用于检查，不参与选择。
# XXX 将sample和holdout分别写成只含目标field的/tmp NekRS临时文件，不会在原数据目录留下文件。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
orders=["XYZ","XZY","YXZ","YZX","ZXY","ZYX"]
datasets=[
    {"name":"test0","path":os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),"field":"T","max_trees":4,"windows":[(0,512,128)]},
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"field":"UZ","max_trees":5,"windows":[(0,1600,400)]},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"field":"UZ","max_trees":5,"windows":[(0,512,128)]},
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","field":"UZ","max_trees":5,"windows":[(1335789,1353197,1024),(1840664,1858072,1024)]},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","field":"UZ","max_trees":5,"windows":[(578353,595761,1024),(1229283,1246691,1024)]},
]
relative_eb=1e-4
detail_path=os.path.join(experiment_dir,"e26_sz3_sampled_configuration_details.csv")
summary_path=os.path.join(experiment_dir,"e26_sz3_sampled_configuration_summary.csv")

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

def cleanup_method_outputs(data_path,field,method):
    remove_if_exists(f"{data_path}_{field}.{method}")
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

detail_fields=["Dataset","Window","Role","Partition Begin","Partition End","Element Count","Abs EB","Order","Predictor","Trees","Method","BR","PSNR","Max Error","Compression Throughput","Decompression Throughput"]
summary_fields=["Dataset","Window","Partition Begin","Partition End","Sample Elements","Holdout Elements","Abs EB","Sample Order","Holdout Order","Sample Predictor","Holdout Predictor","Sample Trees","Holdout Trees","Sample Configuration","Holdout Configuration","Baseline Holdout BR","Sample-selected Holdout BR","Holdout-selected BR","Holdout Gain vs Baseline (%)","Holdout Regret (%)","Sample-selected Holdout PSNR","Holdout-selected PSNR"]
temp_dir=tempfile.mkdtemp(prefix="mrcz-e26-")
summary_rows=[]
try:
    with open(detail_path,"w",newline="") as detail_file:
        detail_writer=csv.DictWriter(detail_file,fieldnames=detail_fields)
        detail_writer.writeheader()
        for dataset in datasets:
            header,original_header,field_values=open_field(dataset["path"],dataset["field"])
            for window_id,(partition_begin,partition_end,sample_size) in enumerate(dataset["windows"],1):
                sample_path=os.path.join(temp_dir,f"{dataset['name']}-{window_id}-sample.f00001")
                holdout_path=os.path.join(temp_dir,f"{dataset['name']}-{window_id}-holdout.f00001")
                relative_indices=np.arange(partition_end-partition_begin,dtype=np.int64)
                rng=np.random.default_rng(2600+window_id)
                sample_relative_indices=np.sort(rng.choice(relative_indices,size=sample_size,replace=False))
                holdout_relative_indices=np.setdiff1d(relative_indices,sample_relative_indices,assume_unique=True)
                sample_indices=sample_relative_indices+partition_begin
                holdout_indices=holdout_relative_indices+partition_begin
                write_compact_nekrs(sample_path,header,original_header,field_values,sample_indices)
                write_compact_nekrs(holdout_path,header,original_header,field_values,holdout_indices)
                role_paths={"sample":sample_path,"holdout":holdout_path}
                field="T"
                partition_min=float(np.min(field_values[partition_begin:partition_end]))
                partition_max=float(np.max(field_values[partition_begin:partition_end]))
                abs_eb=(partition_max-partition_min)*relative_eb
                window_error_args=f"-M ABS {abs_eb:.17g}"
                cache={}
                try:
                    def measure(role,order,predictor,tree_num):
                        key=(role,order,predictor,tree_num)
                        if key in cache:
                            return cache[key]
                        data_path=role_paths[role]
                        method=method_name(order,predictor,tree_num)
                        cleanup_method_outputs(data_path,field,method)
                        base_args=f"-format nekrs -i {data_path} {window_error_args} -method {method} -max-cpus 1 -fields [{field}]"
                        compression_output=run_command(f"./install/bin/mrcz_compress {base_args}",silent=True)
                        decompression_output=run_command(f"./install/bin/mrcz_decompress {base_args}",silent=True)
                        result={
                            "BR":parse_number(decompression_output,r"BR= ([0-9.eE+-]+)","BR"),"PSNR":parse_number(decompression_output,r"PSNR= ([0-9.eE+-]+)","PSNR"),
                            "Max Error":parse_number(decompression_output,r"max_err= ([0-9.eE+-]+)","max error"),"Compression Throughput":parse_number(compression_output,r"compression bandwidth per cpu= ([0-9.eE+-]+)","compression throughput"),
                            "Decompression Throughput":parse_number(decompression_output,r"decompression bandwidth per cpu= ([0-9.eE+-]+)","decompression throughput"),"Method":method,
                        }
                        cache[key]=result
                        detail_writer.writerow({"Dataset":dataset["name"],"Window":window_id,"Role":role,"Partition Begin":partition_begin,"Partition End":partition_end,"Element Count":len(sample_indices) if role=="sample" else len(holdout_indices),"Abs EB":abs_eb,"Order":order,"Predictor":predictor,"Trees":tree_num,**result})
                        detail_file.flush()
                        cleanup_method_outputs(data_path,field,method)
                        return result

                    def select_configuration(role):
                        order=min(orders,key=lambda candidate:measure(role,candidate,"linear",2)["BR"])
                        predictor=min(["linear","cubic"],key=lambda candidate:measure(role,order,candidate,2)["BR"])
                        tree_num=min(range(1,dataset["max_trees"]+1),key=lambda candidate:measure(role,order,predictor,candidate)["BR"])
                        return order,predictor,tree_num

                    sample_configuration=select_configuration("sample")
                    holdout_configuration=select_configuration("holdout")
                    sample_selected_holdout=measure("holdout",*sample_configuration)
                    holdout_selected=measure("holdout",*holdout_configuration)
                    baseline=measure("holdout","XYZ","linear",2)
                    summary_row={
                        "Dataset":dataset["name"],"Window":window_id,"Partition Begin":partition_begin,"Partition End":partition_end,"Sample Elements":len(sample_indices),"Holdout Elements":len(holdout_indices),"Abs EB":abs_eb,
                        "Sample Order":sample_configuration[0],"Holdout Order":holdout_configuration[0],"Sample Predictor":sample_configuration[1],"Holdout Predictor":holdout_configuration[1],
                        "Sample Trees":sample_configuration[2],"Holdout Trees":holdout_configuration[2],"Sample Configuration":"/".join(map(str,sample_configuration)),"Holdout Configuration":"/".join(map(str,holdout_configuration)),
                        "Baseline Holdout BR":baseline["BR"],"Sample-selected Holdout BR":sample_selected_holdout["BR"],"Holdout-selected BR":holdout_selected["BR"],
                        "Holdout Gain vs Baseline (%)":(sample_selected_holdout["BR"]/baseline["BR"]-1)*100,"Holdout Regret (%)":(sample_selected_holdout["BR"]/holdout_selected["BR"]-1)*100,
                        "Sample-selected Holdout PSNR":sample_selected_holdout["PSNR"],"Holdout-selected PSNR":holdout_selected["PSNR"],
                    }
                    summary_rows.append(summary_row)
                    print(f"{dataset['name']} window {window_id}: sample={summary_row['Sample Configuration']},holdout={summary_row['Holdout Configuration']},holdout gain={summary_row['Holdout Gain vs Baseline (%)']:.3f}%,regret={summary_row['Holdout Regret (%)']:.3f}%")
                finally:
                    for role_path in role_paths.values():
                        remove_if_exists(f"{role_path}_header.mrcz")
                        remove_if_exists(f"{role_path}_element_map.mrcz")
    with open(summary_path,"w",newline="") as summary_file:
        writer=csv.DictWriter(summary_file,fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)
finally:
    shutil.rmtree(temp_dir)

print(f"详细结果已写入{detail_path}")
print(f"汇总结果已写入{summary_path}")
