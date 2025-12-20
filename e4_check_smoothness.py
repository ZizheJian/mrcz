import yaml,struct,os,time
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from read_nek5000_header import nek5000_header_c
from py_funcs.quantize_ndarray import quantize_ndarray1,quantize_ndarray2

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

if os.environ.get("SLURM_NTASKS") is not None:
    thread_num=int(os.environ.get("SLURM_NTASKS"))
else:
    thread_num=1

with open(f"{file_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    header=nek5000_header_c()
    header.__dict__.update(header_and_block_info["header"])
    block_num=header_and_block_info["block_num"]
    blocks=header_and_block_info["blocks"]

p=np.zeros((header.element_shape[0],header.element_shape[0]),dtype=np.float64)
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header.element_shape[0]-1).deriv().roots(),[1]))
p[0]=1
p[1]=x
for i in range(2,header.element_shape[0]):
    p[i]=((2*i-1)*x*p[i-1]-(i-1)*p[i-2])/i
w=2/(header.element_shape[0]*(header.element_shape[0]-1)*p[-1]**2)
d=np.einsum("ia,a,ia->i",p,w,p)
element_base_weight=np.ones((8,*header.element_shape),dtype=np.float64)
element_base_weight[0]=np.einsum("a,b,c->abc",1-x,1-x,1-x)/8
element_base_weight[1]=np.einsum("a,b,c->abc",1-x,1-x,x+1)/8
element_base_weight[2]=np.einsum("a,b,c->abc",1-x,x+1,1-x)/8
element_base_weight[3]=np.einsum("a,b,c->abc",1-x,x+1,x+1)/8
element_base_weight[4]=np.einsum("a,b,c->abc",x+1,1-x,1-x)/8
element_base_weight[5]=np.einsum("a,b,c->abc",x+1,1-x,x+1)/8
element_base_weight[6]=np.einsum("a,b,c->abc",x+1,x+1,1-x)/8
element_base_weight[7]=np.einsum("a,b,c->abc",x+1,x+1,x+1)/8

def calc_smoothness(chunk):
    l1_err=np.zeros(header.element_shape[0]-2,dtype=np.float64)
    l2_err=np.zeros(header.element_shape[0]-2,dtype=np.float64)
    for element_id in range(chunk.shape[1]):
        element=chunk[0,element_id]
        err0=element[1:-1,:,:]-(element[:-2,:,:]+element[2:,:,:])/2
        err1=(element[:,1:-1,:]-(element[:,:-2,:]+element[:,2:,:])/2).transpose(1,0,2)
        err2=(element[:,:,1:-1]-(element[:,:,:-2]+element[:,:,2:])/2).transpose(2,0,1)
        l1_err+=np.sum(np.abs(err0)+np.abs(err1)+np.abs(err2),axis=(1,2))
        l2_err+=np.sum(err0**2+err1**2+err2**2,axis=(1,2))
    return l1_err,l2_err

def calc_smoothness_v2(chunk):
    l1_err=np.zeros(header.element_shape[0]-2,dtype=np.float64)
    l2_err=np.zeros(header.element_shape[0]-2,dtype=np.float64)
    for element_id in range(chunk.shape[1]):
        element=chunk[0,element_id]
        err0=element[1:-1,:,:]-(element[:-2,:,:]*(x[2:]-x[1:-1])[:,None,None]+(element[2:,:,:]* (x[1:-1]-x[:-2])[:,None,None]))/(x[2:]-x[:-2])[:,None,None]
        err1=(element[:,1:-1,:]-(element[:,:-2,:]*(x[2:]-x[1:-1])[None,:,None]+(element[:,2:,:]* (x[1:-1]-x[:-2])[None,:,None]))/(x[2:]-x[:-2])[None,:,None]).transpose(1,0,2)
        err2=(element[:,:,1:-1]-(element[:,:,:-2]*(x[2:]-x[1:-1])[None,None,:]+(element[:,:,2:]* (x[1:-1]-x[:-2])[None,None,:]))/(x[2:]-x[:-2])[None,None,:]).transpose(2,0,1)
        l1_err+=np.sum(np.abs(err0)+np.abs(err1)+np.abs(err2),axis=(1,2))
        l2_err+=np.sum(err0**2+err1**2+err2**2,axis=(1,2))
    return l1_err,l2_err

with ProcessPoolExecutor(max_workers=thread_num) as exe:
    l1_err=np.zeros(header.element_shape[0]-2,dtype=np.float64)
    l2_err=np.zeros(header.element_shape[0]-2,dtype=np.float64)
    for block in blocks:
        block_element_num=os.path.getsize(block["pressure_path"])//(header.word_size*header.element_size)
        with open(block["pressure_path"],"rb") as f_pressure:
            cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=block_element_num*header.element_size).reshape(1,-1,*header.element_shape)
        start_element_id=[(block_element_num*i)//thread_num for i in range(thread_num)]
        end_element_id=[(block_element_num*(i+1))//thread_num for i in range(thread_num)]
        exe_input=[cache[:,start_element_id[i]:end_element_id[i]] for i in range(thread_num)]
        results=list(exe.map(calc_smoothness_v2,exe_input))
        l1_err+=np.sum([res[0] for res in results],axis=0)
        l2_err+=np.sum([res[1] for res in results],axis=0)
    l1_err/=3*header.file_element_num*(header.element_shape[0]**2)
    l2_err/=3*header.file_element_num*(header.element_shape[0]**2)
    l2_err=np.sqrt(l2_err)
    print("L1 err=",l1_err)
    print("L2 err=",l2_err)
    print("prediction_distance=",x[2:]-x[:-2])