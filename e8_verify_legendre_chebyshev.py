# TODO: 3阶样条插值、改变预测顺序
# TODO: 对element按照连续性进行分类
# TODO: tune各种参数
# TODO: 所有huffman共用一个字典
# TODO: 根据error bound选择数据类型

import yaml,struct,os,time
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from sklearn.preprocessing import PolynomialFeatures
from sklearn.svm import SVC
from py_funcs.read_nekrs_header import nekrs_header_c
from py_funcs.quantize import quantize_ndarray1,quantize_ndarray2,quantize_ndarray3

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

if os.environ.get("SLURM_NTASKS") is not None:
    thread_num=int(os.environ.get("SLURM_NTASKS"))
else:
    thread_num=1

with open(f"{file_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    header=nekrs_header_c()
    header.__dict__.update(header_and_block_info["header"])
    data_range=header_and_block_info["data_range"]
    block_num=header_and_block_info["block_num"]
    # block_num=1
    blocks=header_and_block_info["blocks"]
    # blocks=header_and_block_info["blocks"][0:block_num]

x_legendre=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header.element_shape[0]-1).deriv().roots(),[1]))
x_chebyshev=np.sort(np.cos((np.arange(header.element_shape[0]))*np.pi/(header.element_shape[0]-1)))
print("Legendre points:",x_legendre)
print("Chebyshev points:",x_chebyshev)

def sum_point_position(chunk):
    point_position=np.zeros((header.element_shape[0],),dtype=np.float64)
    for element_id in range(chunk.shape[1]):
        element=chunk[:,element_id]
        temp0=((element[:,1:,:,:]-element[:,:-1,:,:])**2).sum(axis=0)
        temp0=2*temp0/temp0.sum(axis=0,keepdims=True)
        temp0=np.concatenate((np.zeros((1,header.element_shape[1],header.element_shape[2])),temp0),axis=0)
        temp0=np.cumsum(temp0,axis=0)-1
        temp1=((element[:,:,1:,:]-element[:,:,:-1,:])**2).sum(axis=0)
        temp1=2*temp1/temp1.sum(axis=1,keepdims=True)
        temp1=np.concatenate((np.zeros((header.element_shape[0],1,header.element_shape[2])),temp1),axis=1)
        temp1=np.cumsum(temp1,axis=1)-1
        temp2=((element[:,:,:,1:]-element[:,:,:,:-1])**2).sum(axis=0)
        temp2=2*temp2/temp2.sum(axis=2,keepdims=True)
        temp2=np.concatenate((np.zeros((header.element_shape[0],header.element_shape[1],1)),temp2),axis=2)
        temp2=np.cumsum(temp2,axis=2)-1
        point_position+=np.sum(temp0,axis=(1,2))+np.sum(temp1,axis=(0,2))+np.sum(temp2,axis=(0,1))
    return point_position


with ProcessPoolExecutor(max_workers=thread_num) as exe:
    for block in blocks:
        point_position=np.zeros((header.element_shape[0],),dtype=np.float64)
        print("Processing block:",block["block_id"])
        block_element_num=block["element_num"]
        with open(block["coordinate_path"],"rb") as f_coordinate:
            cache=np.fromfile(f_coordinate,dtype=header.np_dtype,count=3*block_element_num*header.element_size).reshape(3,-1,*header.element_shape).astype(np.float64)
        start_element_id=[(block_element_num*i)//thread_num for i in range(thread_num)]
        end_element_id=[(block_element_num*(i+1))//thread_num for i in range(thread_num)]
        exe_input=[cache[:,start_element_id[i]:end_element_id[i]] for i in range(thread_num)]
        results=list(exe.map(sum_point_position,exe_input))
        for local_point_position in results:
            point_position+=local_point_position

        point_position/=(block_element_num*(header.element_shape[0]**2)*3)
        print("Point position:",point_position)
        exit()
