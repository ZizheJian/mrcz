import yaml
import numpy as np

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
block_id=0
element_id=0

with open(f"{file_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    header=header_and_block_info["header"]
    block_num=header_and_block_info["block_num"]
    blocks=header_and_block_info["blocks"]
block=blocks[block_id]

with open(block["pressure_path"],"rb") as f:
    pressure=np.fromfile(f,dtype=header["np_dtype"]).reshape((1,-1)+header["element_shape"])
    print(f"pressure shape: {pressure.shape}")
pressure=pressure[:,element_id:element_id+1]
pressure=pressure.astype(np.float64)
center=pressure[0,0,4:6,4:6,4:6].sum()/8
pressure-=center
print(f"pressure max={pressure.max()}, min={pressure.min()}")

Ch=np.zeros((header["element_shape"][0],header["element_shape"][0]),dtype=pressure.dtype)
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header["element_shape"][0]-1).deriv().roots(),[1]))
Ch[0]=1
Ch[1]=x
for i in range(2,header["element_shape"][0]):
    Ch[i]=2*x*Ch[i-1]-Ch[i-2]
neg_Ch=np.linalg.inv(Ch)

f=np.einsum("abc,ai->ibc",pressure[0,0],neg_Ch)
f=np.einsum("ibc,bj->ijc",f,neg_Ch)
f=np.einsum("ijc,ck->ijk",f,neg_Ch)
print(f"f.max()={f.max()}, f.min()={f.min()}")

pressure2=np.einsum("ijk,kc->ijc",f,Ch)
pressure2=np.einsum("ijc,jb->ibc",pressure2,Ch)
pressure2=np.einsum("ibc,ia->abc",pressure2,Ch)

print(f"max_err= {np.abs(pressure2-pressure[0,0]).max()}")