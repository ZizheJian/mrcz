# 功能：比较同element角点Lorenzo预测与前一个编号element同位置角点预测。
# 用法：直接运行python3 e10_compare_corner_predictability.py。默认测试5个NekRS数据的UZ，REL error bound为1e-4；大数据各抽取3段4096个element。
# 同element方案按000、001、010、011、100、101、110、111顺序，只使用已经重建的角点；跨element方案用前一个编号element的同位置重建角点。

import math,os
from collections import Counter
import numpy as np
from py_funcs.read_nekrs_header import read_nekrs_header

rel_eb=1e-4
experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
datasets=(
    ("test0",os.path.join(project_root,"third_party/fld_pyutils/demos/data/test0.f00001"),None,None),
    ("mixlay0",os.path.join(project_root,"third_party/pySEMTools/examples/data/mixlay0.f00001"),None,None),
    ("tc_channel0",os.path.join(project_root,"third_party/pySEMTools/examples/data/tc_channel0.f00001"),None,None),
    ("pb0","/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001",((666618,670714),(955020,959116),(1398182,1402278)),16.83629458334),
    ("jicf0","/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000",((294168,298264),(1139024,1143120),(2674406,2678502)),10.4476763947),
)

def round_away_from_zero(values):
    return np.copysign(np.floor(np.abs(values)+0.5),values).astype(np.int64)

def open_uz(path):
    with open(path,"rb") as f:
        header=read_nekrs_header(f)
    offset=136+header.file_element_num*4
    value_num=header.file_element_num*header.element_size
    for variable in header.variables:
        if variable=="X":
            offset+=3*value_num*header.word_size
        elif variable=="U":
            dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
            data=np.memmap(path,dtype=dtype,mode="r",offset=offset,shape=(header.file_element_num,3,header.element_size))
            return header,data[:,2,:]
        elif variable in ("P","T"):
            offset+=value_num*header.word_size
    raise ValueError(f"{path}不包含速度场U")

def corner_indices(shape):
    n0,n1,n2=shape
    return np.array([(i0*n1+i1)*n2+i2 for i0 in ((0,n0-1) if n0>1 else (0,)) for i1 in ((0,n1-1) if n1>1 else (0,)) for i2 in ((0,n2-1) if n2>1 else (0,))],dtype=np.int64)

def lorenzo_prediction(values,index,corner_shape):
    coordinates=np.unravel_index(index,corner_shape)
    active_axes=[axis for axis,coordinate in enumerate(coordinates) if coordinate==1]
    prediction=0.0
    for mask in range(1,1<<len(active_axes)):
        neighbor=list(coordinates)
        selected_num=0
        for bit,axis in enumerate(active_axes):
            if mask&(1<<bit):
                neighbor[axis]=0
                selected_num+=1
        sign=1 if selected_num%2==1 else -1
        prediction+=sign*values[np.ravel_multi_index(tuple(neighbor),corner_shape)]
    return prediction

class ErrorStats:
    def __init__(self):
        self.count=0
        self.abs_sum=0.0
        self.square_sum=0.0

    def add(self,errors):
        errors=np.asarray(errors,dtype=np.float64)
        self.count+=errors.size
        self.abs_sum+=np.abs(errors).sum(dtype=np.float64)
        self.square_sum+=np.square(errors).sum(dtype=np.float64)

    def result(self,data_range):
        return self.abs_sum/self.count/data_range,math.sqrt(self.square_sum/self.count)/data_range

class QbStats:
    def __init__(self):
        self.count=0
        self.zero_count=0
        self.frequencies=Counter()

    def add(self,qb):
        values,counts=np.unique(qb,return_counts=True)
        self.count+=qb.size
        self.zero_count+=int(np.count_nonzero(qb==0))
        self.frequencies.update({int(value):int(count) for value,count in zip(values,counts)})

    def result(self):
        entropy=0.0
        for count in self.frequencies.values():
            probability=count/self.count
            entropy-=probability*math.log2(probability)
        return self.zero_count/self.count,entropy

def evaluate_dataset(name,path,ranges,data_range):
    header,data=open_uz(path)
    if ranges is None:
        ranges=((0,header.file_element_num),)
    indices=corner_indices(header.element_shape)
    corner_shape=tuple(2 if n>1 else 1 for n in header.element_shape)
    chunks=[np.asarray(data[begin:end,indices],dtype=np.float64) for begin,end in ranges]
    if data_range is None:
        data_range=float(data.max()-data.min())
    step=2*rel_eb*data_range
    within_error=ErrorStats()
    adjacent_error=ErrorStats()
    within_qb_stats=QbStats()
    adjacent_qb_stats=QbStats()
    hybrid_qb_stats=QbStats()

    for truth in chunks:
        within_qb=np.empty(truth.shape,dtype=np.int64)
        for element_id in range(truth.shape[0]):
            reconstructed=np.empty(truth.shape[1],dtype=np.float64)
            for corner_id in range(truth.shape[1]):
                if corner_id==0:
                    prediction=0.0
                else:
                    prediction=lorenzo_prediction(reconstructed,corner_id,corner_shape)
                    original_prediction=lorenzo_prediction(truth[element_id],corner_id,corner_shape)
                    within_error.add(original_prediction-truth[element_id,corner_id])
                qb=round_away_from_zero(np.array([(prediction-truth[element_id,corner_id])/step]))[0]
                within_qb[element_id,corner_id]=qb
                reconstructed[corner_id]=prediction-qb*step
        within_qb_stats.add(within_qb)

        adjacent_error.add(truth[:-1]-truth[1:])
        adjacent_qb=np.empty(truth.shape,dtype=np.int64)
        previous=np.zeros(truth.shape[1],dtype=np.float64)
        for element_id in range(truth.shape[0]):
            prediction=np.zeros(truth.shape[1],dtype=np.float64) if element_id==0 else previous
            qb=round_away_from_zero((prediction-truth[element_id])/step)
            adjacent_qb[element_id]=qb
            previous=prediction-qb*step
        adjacent_qb_stats.add(adjacent_qb)

        hybrid_qb=np.empty(truth.shape,dtype=np.int64)
        previous_anchor=0.0
        for element_id in range(truth.shape[0]):
            reconstructed=np.empty(truth.shape[1],dtype=np.float64)
            for corner_id in range(truth.shape[1]):
                prediction=previous_anchor if corner_id==0 else lorenzo_prediction(reconstructed,corner_id,corner_shape)
                qb=round_away_from_zero(np.array([(prediction-truth[element_id,corner_id])/step]))[0]
                hybrid_qb[element_id,corner_id]=qb
                reconstructed[corner_id]=prediction-qb*step
            previous_anchor=reconstructed[0]
        hybrid_qb_stats.add(hybrid_qb)

    within_nmae,within_nrmse=within_error.result(data_range)
    adjacent_nmae,adjacent_nrmse=adjacent_error.result(data_range)
    within_zero,within_entropy=within_qb_stats.result()
    adjacent_zero,adjacent_entropy=adjacent_qb_stats.result()
    hybrid_zero,hybrid_entropy=hybrid_qb_stats.result()
    sampled_elements=sum(end-begin for begin,end in ranges)
    print(f"{name}: shape={header.element_shape}, sampled_elements={sampled_elements}, global_UZ_range={data_range:.9e}")
    print(f"  within-element Lorenzo: NMAE={within_nmae:.9e}, NRMSE={within_nrmse:.9e}, zero_qb={within_zero*100:.3f}%, entropy={within_entropy:.6f} bits/corner")
    print(f"  previous-element same corner: NMAE={adjacent_nmae:.9e}, NRMSE={adjacent_nrmse:.9e}, zero_qb={adjacent_zero*100:.3f}%, entropy={adjacent_entropy:.6f} bits/corner")
    print(f"  hybrid previous anchor + within Lorenzo: zero_qb={hybrid_zero*100:.3f}%, entropy={hybrid_entropy:.6f} bits/corner")
    print(f"  adjacent/within: NMAE={adjacent_nmae/within_nmae:.6f}, NRMSE={adjacent_nrmse/within_nrmse:.6f}, entropy={adjacent_entropy/within_entropy:.6f}")

for dataset in datasets:
    evaluate_dataset(*dataset)
