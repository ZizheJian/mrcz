import csv,os,re,shutil,tempfile,time
import numpy as np
from py_funcs.read_nekrs_header import read_nekrs_header
from third_party.run_command.run_command import run_command

# XXX 比较固定SZ3配置与按partition、按element局部位置采样得到的residual bias。
# XXX bias量化单位为abs_eb/8，每个局部位置搜索宽度为2*abs_eb的最密窗口，并使用最小密度和改善幅度阈值；角点也单独学习。
# XXX bias表写入每个partition的Zstd流，因此报告的BR已经包含metadata。
# XXX 五个case各代表一个NekRS数据；大数据只复制一个连续partition窗口到/tmp，原数据目录不会产生文件。

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
cases=[
    {"name":"test0","path":os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),"field":"T","rel_eb":1e-4,"order":"XYZ","predictor":"linear","trees":4,"window":(0,512)},
    {"name":"mixlay0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),"field":"T","rel_eb":1e-4,"order":"XYZ","predictor":"cubic","trees":5,"window":(0,1600)},
    {"name":"tc_channel0","path":os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),"field":"UZ","rel_eb":1e-4,"order":"YZX","predictor":"cubic","trees":5,"window":(0,512)},
    {"name":"pb0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001","field":"UZ","rel_eb":1e-4,"order":"ZYX","predictor":"cubic","trees":5,"window":(1335789,1353197)},
    {"name":"jicf0","path":"/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000","field":"UZ","rel_eb":1e-6,"order":"ZXY","predictor":"cubic","trees":5,"window":(578353,595761)},
]
result_path=os.path.join(experiment_dir,"e28_sz3_position_bias_results.csv")

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

def write_compact_nekrs(path,header,original_header,values,begin,end):
    element_num=end-begin
    header_terms=original_header[:132].split()
    header_terms[5]=str(element_num).encode()
    header_terms[6]=str(element_num).encode()
    header_terms[11]=b"T"
    new_header=b" ".join(header_terms)
    if len(new_header)>132:
        raise RuntimeError("重写后的NekRS header超过132 bytes")
    map_dtype=np.dtype(np.int32).newbyteorder(header.endian_mode)
    with open(path,"wb") as output_file:
        output_file.write(new_header.ljust(132,b" "))
        output_file.write(original_header[132:136])
        np.arange(1,element_num+1,dtype=map_dtype).tofile(output_file)
        np.asarray(values[begin:end]).tofile(output_file)

def measure(data_path,abs_eb,method):
    cleanup_outputs(data_path,method)
    args=f"-format nekrs -i {data_path} -M ABS {abs_eb:.17g} -method {method} -max-cpus 1 -fields [T]"
    compression_output=run_command(f"./install/bin/mrcz_compress {args}",silent=True)
    decompression_output=run_command(f"./install/bin/mrcz_decompress {args}",silent=True)
    result={
        "BR":parse_number(decompression_output,r"BR= ([0-9.eE+-]+)","BR"),"PSNR":parse_number(decompression_output,r"PSNR= ([0-9.eE+-]+)","PSNR"),
        "Max Error":parse_number(decompression_output,r"max_err= ([0-9.eE+-]+)","max error"),"Compression Throughput":parse_number(compression_output,r"compression bandwidth per cpu= ([0-9.eE+-]+)","compression throughput"),
        "Decompression Throughput":parse_number(decompression_output,r"decompression bandwidth per cpu= ([0-9.eE+-]+)","decompression throughput"),
    }
    cleanup_outputs(data_path,method)
    return result

rows=[]
temp_dir=tempfile.mkdtemp(prefix="mrcz-e28-")
experiment_start=time.perf_counter()
try:
    for case in cases:
        header,original_header,values=open_field(case["path"],case["field"])
        begin,end=case["window"]
        data_path=os.path.join(temp_dir,f"{case['name']}.f00001")
        write_compact_nekrs(data_path,header,original_header,values,begin,end)
        data_range=float(np.max(values[begin:end])-np.min(values[begin:end]))
        abs_eb=data_range*case["rel_eb"]
        suffix=f"{case['order']}-{case['predictor']}-{case['trees']}"
        baseline_method=f"SZ3-3D-GLL-sampled-{suffix}"
        bias_method=f"SZ3-3D-GLL-position-bias-{suffix}"
        baseline=measure(data_path,abs_eb,baseline_method)
        bias=measure(data_path,abs_eb,bias_method)
        row={"Dataset":case["name"],"Field":case["field"],"Elements":end-begin,"REL EB":case["rel_eb"],"Abs EB":abs_eb,"Order":case["order"],"Predictor":case["predictor"],"Trees":case["trees"],
            "Baseline BR":baseline["BR"],"Bias BR":bias["BR"],"BR Change (%)":(bias["BR"]/baseline["BR"]-1)*100,"Baseline PSNR":baseline["PSNR"],"Bias PSNR":bias["PSNR"],
            "Baseline Max Error":baseline["Max Error"],"Bias Max Error":bias["Max Error"],"Baseline Max Error / EB":baseline["Max Error"]/abs_eb,"Bias Max Error / EB":bias["Max Error"]/abs_eb,"Baseline Compression Throughput":baseline["Compression Throughput"],"Bias Compression Throughput":bias["Compression Throughput"],
            "Baseline Decompression Throughput":baseline["Decompression Throughput"],"Bias Decompression Throughput":bias["Decompression Throughput"]}
        rows.append(row)
        print(f"{case['name']} {case['field']}: BR {baseline['BR']:.6f}->{bias['BR']:.6f} ({row['BR Change (%)']:+.3f}%),PSNR {baseline['PSNR']:.3f}->{bias['PSNR']:.3f}",flush=True)
        remove_if_exists(f"{data_path}_header.mrcz")
        remove_if_exists(f"{data_path}_element_map.mrcz")
        del values
finally:
    shutil.rmtree(temp_dir)

with open(result_path,"w",newline="") as output_file:
    writer=csv.DictWriter(output_file,fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"实验用时{time.perf_counter()-experiment_start:.3f}秒，结果已写入{result_path}")
