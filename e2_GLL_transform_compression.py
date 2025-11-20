import yaml,struct,psutil,os
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from read_nek5000_header import nek5000_header_c

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
rel_eb=1e-3

if os.environ.get("SLURM_NTASKS") is not None:
    thread_num=int(os.environ.get("SLURM_NTASKS"))
else:
    thread_num=1

def quantize_ndarray(parameter:np.ndarray,eb:float)->np.ndarray:
    mask_positive=(parameter>=0)
    mask_negative=(~mask_positive)
    parameter[mask_positive]=np.ceil((parameter[mask_positive]-eb)/(2*eb))
    parameter[mask_negative]=np.floor((parameter[mask_negative]+eb)/(2*eb))
    return parameter.astype(np.int32),2*eb*parameter

with open(f"{file_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    header=nek5000_header_c()
    header.__dict__.update(header_and_block_info["header"])
    block_num=header_and_block_info["block_num"]
    blocks=header_and_block_info["blocks"]

abs_eb=(header.data_range["P"][1]-header.data_range["P"][0])*rel_eb
param_eb=2*rel_eb*0.5
max_err=0
mse=0
print(f"presure_min={header.data_range['P'][0]}, presure_max={header.data_range['P'][1]}")
print(f"abs_eb={abs_eb}")
print(f"param_eb={param_eb}")

p=np.zeros((header.element_shape[0],header.element_shape[0]),dtype=np.float64)
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header.element_shape[0]-1).deriv().roots(),[1]))
p[0]=1
p[1]=x
for i in range(2,header.element_shape[0]):
    p[i]=((2*i-1)*x*p[i-1]-(i-1)*p[i-2])/i
w=2/(header.element_shape[0]*(header.element_shape[0]-1)*p[-1]**2)
d=np.einsum("ia,a,ia->i",p,w,p)

f_qb_list=[]
qb_list=[]
for block in blocks:
    with open(block["pressure_path"],"rb") as f:
        data=np.fromfile(f,dtype=header.np_dtype).reshape((1,-1)+header.element_shape)
        print(f"data shape: {data.shape}")
    for element_id in range(data.shape[1]):
        print(f"block {block['block_id']+1}/{block_num}, element {element_id+1}/{data.shape[1]}")
        element=data[0,element_id].astype(np.float64)
        if header.element_shape[0]%2==0:
            center=element[header.element_shape[0]//2-1:header.element_shape[0]//2+1,
                           header.element_shape[1]//2-1:header.element_shape[1]//2+1,
                           header.element_shape[2]//2-1:header.element_shape[2]//2+1].sum()/8
        else:
            center=element[header.element_shape[0]//2,
                           header.element_shape[1]//2,
                           header.element_shape[2]//2]
        element-=center

        f=np.einsum("abc,ia,a->ibc",element,p,w)/d.reshape(-1,1,1)
        f=np.einsum("ibc,jb,b->ijc",f,p,w)/d.reshape(1,-1,1)
        f=np.einsum("ijc,kc,c->ijk",f,p,w)/d.reshape(1,1,-1)
        print(f"f.max()={f.max()}, f.min()={f.min()}")

        f_qb,f_quantized=quantize_ndarray(f.copy(),param_eb)
        f_qb_list.append(f_qb)

        pred=np.einsum("ijk,kc->ijc",f_quantized,p)
        pred=np.einsum("ijc,jb->ibc",pred,p)
        pred=np.einsum("ibc,ia->abc",pred,p)
        
        qb,new_element=quantize_ndarray(element-pred,abs_eb)
        qb_list.append(qb)
        new_element+=pred
        
        max_err=max(max_err,np.abs(new_element-element).max())
        mse+=np.sum((new_element-element)**2)

######## Huffman f ########
f_hf_bs=b""
f_values,f_counts=np.unique(f_qb,return_counts=True)
f_freq=dict(zip(f_values,f_counts))
f_hf_bs+=struct.pack("<i",len(f_values))
max_f=np.abs(f_values).max().item()
if max_f<=2**15-1:
    f_hf_bs+=struct.pack("<b",15)
    f_hf_bs+=struct.pack(f"<{len(f_values)}h",*f_values.tolist())
else:
    f_hf_bs+=struct.pack("<b",31)
    f_hf_bs+=struct.pack(f"<{len(f_values)}i",*f_values.tolist())
f_hf_bs+=struct.pack(f"<{len(f_values)}i",*f_counts.tolist())
f_codec=HuffmanCodec.from_frequencies(f_freq)
temp_bs=f_codec.encode(f_qb.flatten().tolist())
f_hf_bs+=struct.pack("<i",len(temp_bs))
f_hf_bs+=temp_bs
print(f"频率编码长度: {len(f_hf_bs)}")
######## Huffman qb ########
hf_bs=b""
qb_values,qb_counts=np.unique(qb,return_counts=True)
qb_freq=dict(zip(qb_values,qb_counts))
hf_bs+=struct.pack("<i",len(qb_values))
max_qb=np.abs(qb_values).max().item()
if max_qb<=2**15-1:
    hf_bs+=struct.pack("<b",15)
    hf_bs+=struct.pack(f"<{len(qb_values)}h",*qb_values.tolist())
else:
    hf_bs+=struct.pack("<b",31)
    hf_bs+=struct.pack(f"<{len(qb_values)}i",*qb_values.tolist())
hf_bs+=struct.pack(f"<{len(qb_values)}i",*qb_counts.tolist())
qb_codec=HuffmanCodec.from_frequencies(qb_freq)
temp_bs=qb_codec.encode(qb.flatten().tolist())
hf_bs+=struct.pack("<i",len(temp_bs))
hf_bs+=temp_bs
print(f"qb编码长度: {len(f_hf_bs)}")
######## Zstd ########
cctx=zstd.ZstdCompressor()
zstd_bs=cctx.compress(f_hf_bs+hf_bs)
print(f"Zstd编码长度: {len(zstd_bs)}")
print(f"原始数据大小: {header.file_element_num*header.element_size*header.word_size} bytes")
print(f"压缩比: {header.file_element_num*header.element_size*header.word_size/len(zstd_bs)}")

print(f"max_err= {max_err}({max_err*100/abs_eb}%%)")
mse/=(header.file_element_num*header.element_size)
psnr=10*np.log10((header.data_range["P"][1]-header.data_range["P"][0])**2/mse)
print(f"MSE={mse}, PSNR={psnr}dB")