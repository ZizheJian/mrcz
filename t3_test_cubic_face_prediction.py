# 功能：在相同的共享面目标点上比较跨element与沿共享面的cubic预测。
# 对每个8x8共享面使用中央4x4个点。跨面预测使用两个element中各两层内部点，
# 并根据到目标的有符号物理距离计算cubic Lagrange权重；沿面预测分别沿两个
# reference-space方向使用左右各两个点，并将两个方向的误差合并统计。

import csv
import os
import re
import numpy as np
from py_funcs.read_nekrs_header import read_nekrs_header

data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001"
requested_fields=("P",)
random_seed=0
sample_element_num=100
neighbor_scan_chunk_size=20000
field_range_chunk_size=20000
face_center_tol=1e-7
result_directory="results/prediction_locality/P"
entropy_relative_eb=1e-4

# 角点顺序为000、001、010、011、100、101、110、111。
face_corner_ids=((0,1,2,3),(0,1,4,5),(0,2,4,6),(1,3,5,7),(2,3,6,7),(4,5,6,7))
# 每个面对应(法向轴,低侧/高侧)。
face_axis_and_side=((0,0),(1,0),(2,0),(2,1),(1,1),(0,1))

def open_nekrs_fields(path,field_names):
    def should_load(field_name):
        return field_names is None or field_name in field_names

    with open(path,"rb") as input_file:
        header=read_nekrs_header(input_file)
        input_file.seek(header.file_element_num*4,os.SEEK_CUR)
        field_offset=input_file.tell()

    dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
    element_num=header.file_element_num
    element_size=header.element_size
    itemsize=dtype.itemsize
    loaded_fields={}
    coordinates=None
    variable_index=0
    while variable_index<len(header.variables):
        variable=header.variables[variable_index]
        if variable=="X":
            raw=np.memmap(path,dtype=dtype,mode="r",offset=field_offset,shape=(element_num,3,element_size))
            coordinates=raw.transpose(0,2,1)
            field_offset+=element_num*3*element_size*itemsize
        elif variable=="U":
            raw=np.memmap(path,dtype=dtype,mode="r",offset=field_offset,shape=(element_num,3,element_size))
            for component_id,name in enumerate(("UX","UY","UZ")):
                if should_load(name):
                    loaded_fields[name]=raw[:,component_id,:]
            field_offset+=element_num*3*element_size*itemsize
        elif variable in ("P","T"):
            if should_load(variable):
                loaded_fields[variable]=np.memmap(path,dtype=dtype,mode="r",offset=field_offset,shape=(element_num,element_size))
            field_offset+=element_num*element_size*itemsize
        elif variable=="S":
            scalar_match=re.match(r"S(\d{2})",header.variables[variable_index:])
            if scalar_match is None:
                raise ValueError(f"无法解析被动标量说明：{header.variables}")
            scalar_num=int(scalar_match.group(1))
            for scalar_id in range(scalar_num):
                name=f"S{scalar_id:02d}"
                if should_load(name):
                    loaded_fields[name]=np.memmap(path,dtype=dtype,mode="r",offset=field_offset,shape=(element_num,element_size))
                field_offset+=element_num*element_size*itemsize
            variable_index+=2
        else:
            raise ValueError(f"不支持的NekRS变量说明：{variable}")
        variable_index+=1

    if coordinates is None:
        raise ValueError("文件不包含坐标场X")
    if not loaded_fields:
        raise ValueError("文件不包含任何指定field")
    if field_names is not None:
        missing_fields=[field_name for field_name in field_names if field_name not in loaded_fields]
        if missing_fields:
            print(f"警告：文件不包含这些指定field，将跳过：{missing_fields}")
    if os.path.getsize(path)<field_offset:
        raise ValueError("文件长度小于header所描述的数据长度")
    print(f"本次测试field：{list(loaded_fields)}")
    return header,coordinates,loaded_fields

def get_corner_positions(n):
    n2=n*n
    return np.array((0,n-1,(n-1)*n,n2-1,(n-1)*n2,(n-1)*n2+n-1,(n-1)*n2+(n-1)*n,n**3-1),dtype=np.int64)

def quantize_face_centers(face_centers):
    return np.rint(face_centers/face_center_tol).astype(np.int64)

def hash_quantized_centers(quantized_centers):
    unsigned=quantized_centers.astype(np.uint64,copy=False)
    return ((unsigned[...,0]*np.uint64(0x9E3779B185EBCA87))^(unsigned[...,1]*np.uint64(0xC2B2AE3D27D4EB4F))^(unsigned[...,2]*np.uint64(0x165667B19E3779F9)))

def find_sample_neighbors(coordinates,corner_positions,sample_element_ids):
    sample_corners=np.asarray(coordinates[sample_element_ids[:,None],corner_positions[None,:],:])
    sample_face_centers=np.stack([np.mean(sample_corners[:,face_ids,:],axis=1,dtype=np.float64) for face_ids in face_corner_ids],axis=1)
    sample_keys=quantize_face_centers(sample_face_centers)
    target_faces={}
    for sample_index,element_id in enumerate(sample_element_ids):
        for face_id in range(6):
            key=tuple(sample_keys[sample_index,face_id])
            target_faces.setdefault(key,[]).append((int(element_id),face_id))
    target_hashes=np.array(list({int(value) for value in hash_quantized_centers(sample_keys).ravel()}),dtype=np.uint64)

    neighbors={int(element_id):[] for element_id in sample_element_ids}
    element_num=len(coordinates)
    for element_start in range(0,element_num,neighbor_scan_chunk_size):
        element_stop=min(element_start+neighbor_scan_chunk_size,element_num)
        chunk_corners=np.asarray(coordinates[element_start:element_stop,corner_positions,:])
        chunk_face_centers=np.stack([np.mean(chunk_corners[:,face_ids,:],axis=1,dtype=np.float64) for face_ids in face_corner_ids],axis=1)
        chunk_keys=quantize_face_centers(chunk_face_centers)
        chunk_hashes=hash_quantized_centers(chunk_keys)
        for local_element_id,face_id in np.argwhere(np.isin(chunk_hashes,target_hashes)):
            element_id=element_start+int(local_element_id)
            key=tuple(chunk_keys[local_element_id,face_id])
            for sample_element_id,sample_face_id in target_faces.get(key,()):
                if element_id!=sample_element_id:
                    neighbors[sample_element_id].append((element_id,sample_face_id,int(face_id)))
        print(f"已扫描{element_stop}/{element_num}个element",flush=True)
    for element_id in neighbors:
        neighbors[element_id]=sorted(set(neighbors[element_id]))
    return neighbors

def get_face_layer(element,face_id,depth):
    axis,side=face_axis_and_side[face_id]
    layer_index=depth if side==0 else element.shape[axis]-1-depth
    return np.take(element,layer_index,axis=axis)

def get_face_point_mapping(coordinate_a,face_a,coordinate_b,face_b):
    face_coordinates_a=get_face_layer(coordinate_a,face_a,0).reshape(-1,3)
    face_coordinates_b=get_face_layer(coordinate_b,face_b,0).reshape(-1,3)
    distance_matrix=np.linalg.norm(face_coordinates_a[:,None,:]-face_coordinates_b[None,:,:],axis=-1)
    mapping=np.argmin(distance_matrix,axis=1)
    distances=distance_matrix[np.arange(len(face_coordinates_a)),mapping]
    element_scale=max(np.ptp(coordinate_a,axis=(0,1,2)).max(),np.ptp(coordinate_b,axis=(0,1,2)).max())
    match_tolerance=max(1e-6*element_scale,1e-10)
    if np.max(distances)>match_tolerance or len(np.unique(mapping))!=len(mapping):
        raise ValueError(f"公共面点没有一一重合：最大距离={np.max(distances)}，容差={match_tolerance}")
    return mapping

def get_reference_nodes(n):
    interior=np.polynomial.legendre.Legendre.basis(n-1).deriv().roots()
    return np.concatenate(([-1.0],interior,[1.0]))

def get_lagrange_weights(reference_positions,target_positions):
    reference_positions=np.asarray(reference_positions,dtype=np.float64)
    if reference_positions.shape[-1]!=4:
        raise ValueError("cubic interpolation必须使用四个参考位置")
    target_positions=np.broadcast_to(np.asarray(target_positions,dtype=np.float64),reference_positions.shape[:-1])
    weights=np.ones(reference_positions.shape,dtype=np.float64)
    for reference_id in range(4):
        for other_id in range(4):
            if other_id==reference_id:
                continue
            denominator=reference_positions[...,reference_id]-reference_positions[...,other_id]
            if np.any(np.abs(denominator)<=np.finfo(np.float64).eps):
                raise ValueError("cubic interpolation包含重合的参考位置")
            weights[...,reference_id]*=(target_positions-reference_positions[...,other_id])/denominator
    return weights

class ErrorAccumulator:
    def __init__(self):
        self.absolute_error_sum=0.0
        self.squared_error_sum=0.0
        self.value_num=0
        self.errors=[]

    def add(self,prediction,target):
        error=np.asarray(prediction,dtype=np.float64)-np.asarray(target,dtype=np.float64)
        self.absolute_error_sum+=float(np.sum(np.abs(error)))
        self.squared_error_sum+=float(np.sum(error*error))
        self.value_num+=error.size
        self.errors.append(error.ravel())

    def mae(self):
        return self.absolute_error_sum/self.value_num

    def rmse(self):
        return np.sqrt(self.squared_error_sum/self.value_num)

    def quantized_counts(self,absolute_eb):
        errors=np.concatenate(self.errors)
        scaled=errors/(2*absolute_eb)
        symbols=(np.sign(scaled)*np.floor(np.abs(scaled)+0.5)).astype(np.int64)
        values,counts=np.unique(symbols,return_counts=True)
        return dict(zip(values,counts))

def empirical_entropy(counts):
    frequencies=np.asarray(list(counts.values()),dtype=np.float64)
    probabilities=frequencies/np.sum(frequencies)
    return float(-np.sum(probabilities*np.log2(probabilities)))

def get_global_field_ranges(fields):
    field_ranges={}
    for field_name,field_data in fields.items():
        field_min=np.inf
        field_max=-np.inf
        for element_start in range(0,len(field_data),field_range_chunk_size):
            element_stop=min(element_start+field_range_chunk_size,len(field_data))
            values=np.asarray(field_data[element_start:element_stop])
            field_min=min(field_min,float(np.nanmin(values)))
            field_max=max(field_max,float(np.nanmax(values)))
        if not np.isfinite(field_min) or not np.isfinite(field_max) or field_max<=field_min:
            raise ValueError(f"field {field_name}没有有效的非零范围")
        field_ranges[field_name]=field_max-field_min
        print(f"{field_name}全局范围：[{field_min:.8e},{field_max:.8e}]，range={field_ranges[field_name]:.8e}",flush=True)
    return field_ranges

def get_cross_face_weights(coordinate_a,face_a,coordinate_b,face_b,mapping,n):
    face_coordinate_a=get_face_layer(coordinate_a,face_a,0).reshape(n,n,3)
    face_coordinate_b=get_face_layer(coordinate_b,face_b,0).reshape(-1,3)[mapping].reshape(n,n,3)
    target_coordinates=0.5*(face_coordinate_a+face_coordinate_b)
    first_a=get_face_layer(coordinate_a,face_a,1)
    second_a=get_face_layer(coordinate_a,face_a,2)
    first_b=get_face_layer(coordinate_b,face_b,1).reshape(-1,3)[mapping].reshape(n,n,3)
    second_b=get_face_layer(coordinate_b,face_b,2).reshape(-1,3)[mapping].reshape(n,n,3)
    distance_first_a=np.linalg.norm(first_a-target_coordinates,axis=-1)
    distance_second_a=np.linalg.norm(second_a-target_coordinates,axis=-1)
    distance_first_b=np.linalg.norm(first_b-target_coordinates,axis=-1)
    distance_second_b=np.linalg.norm(second_b-target_coordinates,axis=-1)
    signed_positions=np.stack((-distance_second_a,-distance_first_a,distance_first_b,distance_second_b),axis=-1)
    return get_lagrange_weights(signed_positions,0.0)

def predict_along_face(face_values,reference_nodes,target_indices):
    target_num=len(target_indices)
    first_predictions=np.empty((target_num,target_num),dtype=np.float64)
    second_predictions=np.empty((target_num,target_num),dtype=np.float64)
    for output_index,target_index in enumerate(target_indices):
        reference_indices=np.array((target_index-2,target_index-1,target_index+1,target_index+2))
        weights=get_lagrange_weights(reference_nodes[reference_indices],reference_nodes[target_index])
        first_predictions[output_index,:]=weights@face_values[reference_indices[:,None],target_indices[None,:]]
        second_predictions[:,output_index]=face_values[target_indices[:,None],reference_indices[None,:]]@weights
    return first_predictions,second_predictions

def evaluate_neighbor_pair(element_a,face_a,element_b,face_b,coordinates,fields,n,target_indices,reference_nodes,method_errors,max_shared_differences):
    coordinate_a=np.asarray(coordinates[element_a]).reshape(n,n,n,3)
    coordinate_b=np.asarray(coordinates[element_b]).reshape(n,n,n,3)
    mapping=get_face_point_mapping(coordinate_a,face_a,coordinate_b,face_b)
    cross_weights=get_cross_face_weights(coordinate_a,face_a,coordinate_b,face_b,mapping,n)
    center=np.ix_(target_indices,target_indices)

    for field_name,field_data in fields.items():
        value_a=np.asarray(field_data[element_a]).reshape(n,n,n)
        value_b=np.asarray(field_data[element_b]).reshape(n,n,n)
        face_value_a=get_face_layer(value_a,face_a,0)
        face_value_b=get_face_layer(value_b,face_b,0).reshape(-1)[mapping].reshape(n,n)
        max_shared_differences[field_name]=max(max_shared_differences[field_name],float(np.max(np.abs(face_value_a-face_value_b))))
        face_values=0.5*(face_value_a+face_value_b)
        targets=face_values[center]

        cross_references=np.stack((get_face_layer(value_a,face_a,2),get_face_layer(value_a,face_a,1),get_face_layer(value_b,face_b,1).reshape(-1)[mapping].reshape(n,n),get_face_layer(value_b,face_b,2).reshape(-1)[mapping].reshape(n,n)),axis=-1)
        cross_predictions=np.sum(cross_weights*cross_references,axis=-1)
        method_errors["cross"][field_name].add(cross_predictions[center],targets)

        first_predictions,second_predictions=predict_along_face(face_values,reference_nodes,target_indices)
        method_errors["along"][field_name].add(first_predictions,targets)
        method_errors["along"][field_name].add(second_predictions,targets)

def build_result_rows(method_errors,field_ranges,face_num,target_num_per_face):
    rows=[]
    for method_name,field_errors in method_errors.items():
        normalized_absolute_error_sum=0.0
        normalized_squared_error_sum=0.0
        normalized_value_num=0
        aggregate_counts={}
        for field_name,error in field_errors.items():
            field_range=field_ranges[field_name]
            nmae=error.mae()/field_range
            nrmse=error.rmse()/field_range
            counts=error.quantized_counts(entropy_relative_eb*field_range)
            entropy=empirical_entropy(counts)
            rows.append({"method":method_name,"field":field_name,"face_num":face_num,"target_num_per_face":target_num_per_face,"value_num":error.value_num,"mae":error.mae(),"rmse":error.rmse(),"nmae":nmae,"nrmse":nrmse,"entropy":entropy})
            for symbol,count in counts.items():aggregate_counts[symbol]=aggregate_counts.get(symbol,0)+count
            normalized_absolute_error_sum+=error.absolute_error_sum/field_range
            normalized_squared_error_sum+=error.squared_error_sum/(field_range*field_range)
            normalized_value_num+=error.value_num
        rows.append({"method":method_name,"field":"ALL","face_num":face_num,"target_num_per_face":target_num_per_face,"value_num":normalized_value_num,"mae":"","rmse":"","nmae":normalized_absolute_error_sum/normalized_value_num,"nrmse":np.sqrt(normalized_squared_error_sum/normalized_value_num),"entropy":empirical_entropy(aggregate_counts)})
    return rows

def print_results(rows):
    for row in rows:
        if row["field"]=="ALL":
            print(f'{row["method"]} ALL: NRMSE={row["nrmse"]:.8e},NMAE={row["nmae"]:.8e},entropy={row["entropy"]:.8f} bits/symbol,N={row["value_num"]}')
        else:
            print(f'{row["method"]} {row["field"]}: RMSE={row["rmse"]:.8e},MAE={row["mae"]:.8e},NRMSE={row["nrmse"]:.8e},NMAE={row["nmae"]:.8e},entropy={row["entropy"]:.8f} bits/symbol,N={row["value_num"]}')

def write_results(rows,path):
    os.makedirs(os.path.dirname(path),exist_ok=True)
    field_names=("method","field","face_num","target_num_per_face","value_num","mae","rmse","nmae","nrmse","entropy")
    with open(path,"w",newline="",encoding="utf-8") as output_file:
        writer=csv.DictWriter(output_file,fieldnames=field_names,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"结果已写入{path}")

def run_experiment():
    header,coordinates,fields=open_nekrs_fields(data_path,requested_fields)
    if len(header.c_element_shape)!=3 or len(set(header.c_element_shape))!=1:
        raise ValueError(f"只支持三个方向尺寸相同的3D element：{header.c_element_shape}")
    n=header.c_element_shape[0]
    if n<8:
        raise ValueError(f"中央4x4 cubic预测要求每个方向至少有8个点，当前为{n}")
    target_indices=np.arange(2,n-2,dtype=np.int64)
    if len(target_indices)!=n-4:
        raise RuntimeError("目标索引构造失败")
    reference_nodes=get_reference_nodes(n)
    field_ranges=get_global_field_ranges(fields)
    rng=np.random.default_rng(random_seed)
    sample_element_ids=np.sort(rng.choice(header.file_element_num,size=min(sample_element_num,header.file_element_num),replace=False))
    print(f"element尺寸：{header.c_element_shape}")
    print(f"随机抽取{len(sample_element_ids)}个element")
    neighbors=find_sample_neighbors(coordinates,get_corner_positions(n),sample_element_ids)

    unique_pairs={}
    for element_a,neighbor_list in neighbors.items():
        for element_b,face_a,face_b in neighbor_list:
            pair_key=tuple(sorted((element_a,element_b)))
            unique_pairs.setdefault(pair_key,(element_a,face_a,element_b,face_b))
    if not unique_pairs:
        raise RuntimeError("抽取的element没有找到任何相邻element")
    print(f"共找到{len(unique_pairs)}个唯一共享面")

    method_errors={method_name:{field_name:ErrorAccumulator() for field_name in fields} for method_name in ("cross","along")}
    max_shared_differences={field_name:0.0 for field_name in fields}
    for element_a,face_a,element_b,face_b in unique_pairs.values():
        evaluate_neighbor_pair(element_a,face_a,element_b,face_b,coordinates,fields,n,target_indices,reference_nodes,method_errors,max_shared_differences)

    target_num_per_face=len(target_indices)**2
    rows=build_result_rows(method_errors,field_ranges,len(unique_pairs),target_num_per_face)
    print_results(rows)
    for field_name,difference in max_shared_differences.items():
        print(f"{field_name}共享面两份数值的最大绝对差：{difference:.8e}")
    result_name=f"{os.path.basename(data_path)}_cubic_face_prediction.csv"
    write_results(rows,os.path.join(result_directory,result_name))

if __name__=="__main__":
    run_experiment()
