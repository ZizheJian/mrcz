# 功能：随机抽取NekRS element，寻找共享面的相邻element，并测试两种一阶预测。
# 默认只读取pb0.f00001。测试时可将data_path改为小型test0.f00001。
# 程序使用np.memmap按需读取field，但寻找邻居时仍需顺序扫描所有element的角点。
#
# pb0.f00001实验结果（random_seed=0，随机抽取100个element）：
# - 跨公共面预测：综合NRMSE=1.57282844e-3，综合NMAE=4.94392979e-4。
# - 公共面旁边预测：综合NRMSE=5.11564666e-4，综合NMAE=1.56877498e-4。
# - 后一种方法的综合NRMSE和NMAE分别降低约67.5%和68.3%。UX、UY、UZ、T均明显改善；P的NMAE及p50--p95更低，但NRMSE略高，说明存在少量较大的尾部误差。
# - 两种方法预测的位置不同，因此结论应表述为：在这里测试的具体方案下，跨element预测没有表现出精度优势，但需要额外寻找相邻element。
#
# 论文画图方案：
# - 使用相对绝对误差|prediction-target|/(field全局max-min)的ECDF，不使用直方图。
# - 横轴采用对数坐标，纵轴为误差不超过横轴值的节点比例；曲线越靠左上越好。
# - 使用2x3子图分别显示UX、UY、UZ、P、T，以及一个各field等权的汇总ECDF。
# - 两种方法各画一条曲线；可用竖线标记1e-4、1e-3和1e-2，并在图内简要标注综合NRMSE/NMAE。图注中明确说明两种方法的预测目标位置不同。

import os
import re
import numpy as np
from scipy.spatial import cKDTree
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from py_funcs.read_nekrs_header import read_nekrs_header

# data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/fld_pyutils/demos/data/test0.f00001"
# fields=["UZ","T","S00","S01"]
data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001"
random_seed=0
sample_element_num=100
neighbor_scan_chunk_size=20000
field_range_chunk_size=20000
face_center_tol=1e-7

# 角点顺序为000、001、010、011、100、101、110、111。
face_corner_ids=((0,1,2,3),(0,1,4,5),(0,2,4,6),
                 (1,3,5,7),(2,3,6,7),(4,5,6,7))
# 每个面对应(法向轴, 低侧/高侧)。公共面两侧的element不要求使用相同的面编号。
face_axis_and_side=((0,0),(1,0),(2,0),(2,1),(1,1),(0,1))

def open_nekrs_fields(path):
    requested_fields=globals().get("fields")
    def should_load(field_name):
        return requested_fields is None or field_name in requested_fields

    with open(path,"rb") as f:
        header=read_nekrs_header(f)
        f.seek(header.file_element_num*4,os.SEEK_CUR)
        field_offset=f.tell()

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
            raw=np.memmap(path,dtype=dtype,mode="r",offset=field_offset,
                          shape=(element_num,3,element_size))
            coordinates=raw.transpose(0,2,1)
            field_offset+=element_num*3*element_size*itemsize
        elif variable=="U":
            raw=np.memmap(path,dtype=dtype,mode="r",offset=field_offset,
                          shape=(element_num,3,element_size))
            for component_id,name in enumerate(("UX","UY","UZ")):
                if should_load(name):
                    loaded_fields[name]=raw[:,component_id,:]
            field_offset+=element_num*3*element_size*itemsize
        elif variable in ("P","T"):
            if should_load(variable):
                loaded_fields[variable]=np.memmap(
                    path,dtype=dtype,mode="r",offset=field_offset,
                    shape=(element_num,element_size))
            field_offset+=element_num*element_size*itemsize
        elif variable=="S":
            scalar_match=re.match(r"S(\d{2})",header.variables[variable_index:])
            if scalar_match is None:
                raise ValueError(f"无法解析被动标量说明：{header.variables}")
            scalar_num=int(scalar_match.group(1))
            for scalar_id in range(scalar_num):
                name=f"S{scalar_id:02d}"
                if should_load(name):
                    loaded_fields[name]=np.memmap(
                        path,dtype=dtype,mode="r",offset=field_offset,
                        shape=(element_num,element_size))
                field_offset+=element_num*element_size*itemsize
            variable_index+=2
        else:
            raise ValueError(f"不支持的NekRS变量说明：{variable}")
        variable_index+=1

    if coordinates is None:
        raise ValueError("文件不包含坐标场X")
    if not loaded_fields:
        if requested_fields is None:
            raise ValueError("文件不包含除坐标X以外的field")
        raise ValueError(f"文件不包含任何指定field：{requested_fields}")
    if requested_fields is not None:
        missing_fields=[field_name for field_name in requested_fields
                        if field_name not in loaded_fields]
        if missing_fields:
            print(f"警告：文件不包含这些指定field，将跳过：{missing_fields}")
    print(f"本次测试field：{list(loaded_fields)}")
    if os.path.getsize(path)<field_offset:
        raise ValueError("文件长度小于header所描述的数据长度")
    return header,coordinates,loaded_fields

def get_corner_positions(n):
    n2=n*n
    return np.array((0,n-1,(n-1)*n,n2-1,
                     (n-1)*n2,(n-1)*n2+n-1,
                     (n-1)*n2+(n-1)*n,n**3-1),dtype=np.int64)

def quantize_face_centers(face_centers):
    return np.rint(face_centers/face_center_tol).astype(np.int64)

def hash_quantized_centers(quantized_centers):
    unsigned=quantized_centers.astype(np.uint64,copy=False)
    return ((unsigned[...,0]*np.uint64(0x9E3779B185EBCA87)) ^
            (unsigned[...,1]*np.uint64(0xC2B2AE3D27D4EB4F)) ^
            (unsigned[...,2]*np.uint64(0x165667B19E3779F9)))

def find_sample_neighbors(coordinates,corner_pos,sample_element_ids):
    sample_corners=np.asarray(coordinates[sample_element_ids[:,None],corner_pos[None,:],:])
    sample_face_centers=np.stack(
        [np.mean(sample_corners[:,face_ids,:],axis=1,dtype=np.float64)
         for face_ids in face_corner_ids],axis=1)
    sample_keys=quantize_face_centers(sample_face_centers)
    target_faces={}
    for sample_index,element_id in enumerate(sample_element_ids):
        for face_id in range(6):
            key=tuple(sample_keys[sample_index,face_id])
            target_faces.setdefault(key,[]).append((int(element_id),face_id))
    target_hashes=np.array(
        list({int(value) for value in hash_quantized_centers(sample_keys).ravel()}),
        dtype=np.uint64)

    neighbors={int(element_id):[] for element_id in sample_element_ids}
    element_num=len(coordinates)
    for element_start in range(0,element_num,neighbor_scan_chunk_size):
        element_stop=min(element_start+neighbor_scan_chunk_size,element_num)
        chunk_corners=np.asarray(coordinates[element_start:element_stop,corner_pos,:])
        chunk_face_centers=np.stack(
            [np.mean(chunk_corners[:,face_ids,:],axis=1,dtype=np.float64)
             for face_ids in face_corner_ids],axis=1)
        chunk_keys=quantize_face_centers(chunk_face_centers)
        chunk_hashes=hash_quantized_centers(chunk_keys)
        for local_element_id,face_id in np.argwhere(np.isin(chunk_hashes,target_hashes)):
            element_id=element_start+int(local_element_id)
            key=tuple(chunk_keys[local_element_id,face_id])
            for sample_element_id,sample_face_id in target_faces.get(key,()):
                if element_id!=sample_element_id:
                    neighbors[sample_element_id].append(
                        (element_id,sample_face_id,int(face_id)))
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
    distances,mapping=cKDTree(face_coordinates_b).query(face_coordinates_a,k=1)
    element_scale=max(np.ptp(coordinate_a,axis=(0,1,2)).max(),
                      np.ptp(coordinate_b,axis=(0,1,2)).max())
    match_tolerance=max(1e-6*element_scale,1e-10)
    if np.max(distances)>match_tolerance or len(np.unique(mapping))!=len(mapping):
        raise ValueError(
            f"公共面点没有一一重合：最大距离={np.max(distances)}，容差={match_tolerance}")
    return mapping

def get_reference_nodes(n):
    if n<3:
        raise ValueError("公共面旁边的预测至少要求每个方向有3个点")
    interior=np.polynomial.legendre.Legendre.basis(n-1).deriv().roots()
    return np.concatenate(([-1.0],interior,[1.0]))

def adjacent_prediction_weights(n):
    reference_nodes=get_reference_nodes(n)
    distance_to_face=abs(reference_nodes[1]-reference_nodes[0])
    distance_to_second=abs(reference_nodes[2]-reference_nodes[1])
    inverse_distances=np.array((1/distance_to_face,1/distance_to_second))
    return inverse_distances/np.sum(inverse_distances)

class ErrorAccumulator:
    def __init__(self):
        self.absolute_error_sum=0.0
        self.squared_error_sum=0.0
        self.value_num=0
        self.target_min=np.inf
        self.target_max=-np.inf
        self.absolute_errors=[]

    def add(self,prediction,target):
        target=np.asarray(target,dtype=np.float64)
        error=np.asarray(prediction,dtype=np.float64)-target
        self.absolute_error_sum+=np.sum(np.abs(error))
        self.squared_error_sum+=np.sum(error*error)
        self.value_num+=error.size
        self.target_min=min(self.target_min,np.min(target))
        self.target_max=max(self.target_max,np.max(target))
        self.absolute_errors.append(np.abs(error).ravel())

    def mae(self):
        return self.absolute_error_sum/self.value_num

    def rmse(self):
        return np.sqrt(self.squared_error_sum/self.value_num)

    def relative_errors(self,dynamic_range):
        if dynamic_range<=0:
            return np.full(self.value_num,np.nan)
        return np.concatenate(self.absolute_errors)/dynamic_range

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
        if not np.isfinite(field_min) or not np.isfinite(field_max):
            raise ValueError(f"field {field_name}不包含有效的有限数值范围")
        field_ranges[field_name]=field_max-field_min
        print(f"{field_name}全局范围：[{field_min:.8e}, {field_max:.8e}]，"
              f"range={field_ranges[field_name]:.8e}",flush=True)
    return field_ranges

def evaluate_neighbor_pair(element_a,face_a,element_b,face_b,
                           coordinates,fields,n,method_errors):
    coordinate_a=np.asarray(coordinates[element_a]).reshape(n,n,n,3)
    coordinate_b=np.asarray(coordinates[element_b]).reshape(n,n,n,3)
    mapping_a_to_b=get_face_point_mapping(coordinate_a,face_a,coordinate_b,face_b)
    adjacent_weights=adjacent_prediction_weights(n)

    for field_name,field_data in fields.items():
        value_a=np.asarray(field_data[element_a]).reshape(n,n,n)
        value_b=np.asarray(field_data[element_b]).reshape(n,n,n)

        face_value_a=get_face_layer(value_a,face_a,0).reshape(-1)
        face_value_b=get_face_layer(value_b,face_b,0).reshape(-1)[mapping_a_to_b]
        first_inside_a=get_face_layer(value_a,face_a,1).reshape(-1)
        first_inside_b=get_face_layer(value_b,face_b,1).reshape(-1)[mapping_a_to_b]
        cross_prediction=0.5*(first_inside_a+first_inside_b)
        cross_target=0.5*(face_value_a+face_value_b)
        method_errors["cross"][field_name].add(cross_prediction,cross_target)

        second_inside_a=get_face_layer(value_a,face_a,2).reshape(-1)
        adjacent_prediction_a=(adjacent_weights[0]*face_value_a+
                               adjacent_weights[1]*second_inside_a)
        method_errors["adjacent"][field_name].add(adjacent_prediction_a,first_inside_a)

        second_inside_b=get_face_layer(value_b,face_b,2).reshape(-1)[mapping_a_to_b]
        adjacent_prediction_b=(adjacent_weights[0]*face_value_b+
                               adjacent_weights[1]*second_inside_b)
        method_errors["adjacent"][field_name].add(adjacent_prediction_b,first_inside_b)

def print_errors(method_name,field_errors,field_dynamic_ranges):
    combined=ErrorAccumulator()
    normalized_absolute_error_sum=0.0
    normalized_squared_error_sum=0.0
    normalized_value_num=0
    print(method_name)
    for field_name,error in field_errors.items():
        dynamic_range=field_dynamic_ranges[field_name]
        if dynamic_range>0:
            normalized_rmse=error.rmse()/dynamic_range
            normalized_mae=error.mae()/dynamic_range
            normalized_absolute_error_sum+=error.absolute_error_sum/dynamic_range
            normalized_squared_error_sum+=error.squared_error_sum/(dynamic_range*dynamic_range)
            normalized_value_num+=error.value_num
        else:
            normalized_rmse=np.nan
            normalized_mae=np.nan
        relative_errors=error.relative_errors(dynamic_range)
        if dynamic_range>0:
            quantiles=np.quantile(relative_errors,(0.5,0.9,0.95,0.99,1.0))
        else:
            quantiles=np.full(5,np.nan)
        print(f"  {field_name}: RMSE={error.rmse():.8e}, MAE={error.mae():.8e}, "
              f"NRMSE={normalized_rmse:.8e}, NMAE={normalized_mae:.8e}, "
              f"range={dynamic_range:.8e}, N={error.value_num}")
        print(f"    relative |error|: p50={quantiles[0]:.8e}, p90={quantiles[1]:.8e}, "
              f"p95={quantiles[2]:.8e}, p99={quantiles[3]:.8e}, max={quantiles[4]:.8e}")
        combined.absolute_error_sum+=error.absolute_error_sum
        combined.squared_error_sum+=error.squared_error_sum
        combined.value_num+=error.value_num
    if normalized_value_num:
        combined_normalized_rmse=np.sqrt(normalized_squared_error_sum/normalized_value_num)
        combined_normalized_mae=normalized_absolute_error_sum/normalized_value_num
    else:
        combined_normalized_rmse=np.nan
        combined_normalized_mae=np.nan
    print(f"  all raw: RMSE={combined.rmse():.8e}, MAE={combined.mae():.8e}, N={combined.value_num}")
    print(f"  all normalized: NRMSE={combined_normalized_rmse:.8e}, "
          f"NMAE={combined_normalized_mae:.8e}, N={normalized_value_num}")

def write_relative_error_distribution(method_errors,field_dynamic_ranges):
    field_names=list(field_dynamic_ranges)
    fig=make_subplots(rows=len(field_names),cols=1,
                      subplot_titles=[f"{name}: |prediction error| / global range"
                                      for name in field_names])
    method_labels={"cross":"跨公共面","adjacent":"公共面旁边"}
    method_colors={"cross":"#1f77b4","adjacent":"#ff7f0e"}
    for row,field_name in enumerate(field_names,start=1):
        for method_name in ("cross","adjacent"):
            relative_errors=method_errors[method_name][field_name].relative_errors(
                field_dynamic_ranges[field_name])
            fig.add_trace(go.Histogram(
                x=relative_errors,nbinsx=80,histnorm="probability",
                name=method_labels[method_name],legendgroup=method_name,
                marker_color=method_colors[method_name],opacity=0.6,
                showlegend=(row==1)),row=row,col=1)
        fig.update_xaxes(title_text="relative absolute error",row=row,col=1)
        fig.update_yaxes(title_text="probability",row=row,col=1)
    fig.update_layout(barmode="overlay",height=max(350,260*len(field_names)),
                      margin=dict(l=70,r=30,t=60,b=40))
    os.makedirs("fig",exist_ok=True)
    output_path="fig/interpolation_relative_error_distribution.html"
    fig.write_html(output_path)
    print(f"相对误差分布图已写入{output_path}")

header,coordinates,fields=open_nekrs_fields(data_path)
if len(header.c_element_shape)!=3 or len(set(header.c_element_shape))!=1:
    raise ValueError(f"只支持三个方向尺寸相同的3D element：{header.c_element_shape}")
n=header.c_element_shape[0]
corner_pos=get_corner_positions(n)
rng=np.random.default_rng(random_seed)
sample_element_ids=np.sort(rng.choice(header.file_element_num,
                                      size=min(sample_element_num,header.file_element_num),
                                      replace=False))
print(f"element尺寸：{header.c_element_shape}")
field_dynamic_ranges=get_global_field_ranges(fields)
print(f"随机抽取的element：{tuple(map(int,sample_element_ids))}")
neighbors=find_sample_neighbors(coordinates,corner_pos,sample_element_ids)
for element_id in sample_element_ids:
    print(f"element {element_id}找到{len(neighbors[int(element_id)])}个相邻element")

unique_pairs={}
for element_a,neighbor_list in neighbors.items():
    for element_b,face_a,face_b in neighbor_list:
        pair_key=tuple(sorted((element_a,element_b)))
        unique_pairs.setdefault(pair_key,(element_a,face_a,element_b,face_b))
print(f"共找到{len(unique_pairs)}对唯一相邻element")

method_errors={
    method:{field_name:ErrorAccumulator() for field_name in fields}
    for method in ("cross","adjacent")
}
for element_a,face_a,element_b,face_b in unique_pairs.values():
    evaluate_neighbor_pair(element_a,face_a,element_b,face_b,
                           coordinates,fields,n,method_errors)

if not unique_pairs:
    raise RuntimeError("抽取的element没有找到任何相邻element")
print_errors("跨公共面的一阶插值预测",method_errors["cross"],field_dynamic_ranges)
print_errors("公共面旁边的一阶插值预测",method_errors["adjacent"],field_dynamic_ranges)
write_relative_error_distribution(method_errors,field_dynamic_ranges)
