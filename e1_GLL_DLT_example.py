import yaml,struct
import numpy as np
import zstandard as zstd
from plotly import graph_objects as go
from matplotlib import pyplot as plt
from dahuffman import HuffmanCodec

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
block_id=0
element_id=0
lossy=False

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

L=np.zeros((header["element_shape"][0],header["element_shape"][0]),dtype=pressure.dtype)
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header["element_shape"][0]-1).deriv().roots(),[1]))
L[0]=1
L[1]=x
for i in range(2,header["element_shape"][0]):
    L[i]=((2*i-1)*x*L[i-1]-(i-1)*L[i-2])/i
w=2/(header["element_shape"][0]*(header["element_shape"][0]-1)*L[-1]**2)
d=np.einsum("ia,a,ia->i",L,w,L)
# d=2/(2*np.arange(header["element_shape"][0])+1)
print("p=",L)
print("w=",w)
print("d=",d)

f=np.einsum("abc,ia,a->ibc",pressure[0,0],L,w)/d.reshape(-1,1,1)
f=np.einsum("ibc,jb,b->ijc",f,L,w)/d.reshape(1,-1,1)
f=np.einsum("ijc,kc,c->ijk",f,L,w)/d.reshape(1,1,-1)
print(f"f.max()={f.max()}, f.min()={f.min()}")

# fig,axes=plt.subplots(2,5,figsize=(50,20))
# axes=axes.flatten()
# for i in range(10):
#     im=axes[i].imshow(f[i,:,:],origin="lower",cmap="jet",vmin=-1.5,vmax=1.5)
#     axes[i].axis("off")
#     fig.colorbar(im,ax=axes[i])
# plt.savefig(f"temp{element_id}.png")

if lossy:
    def quantize_ndarray(parameter:np.ndarray,eb:float)->np.ndarray:
        mask_positive=(parameter>=0)
        mask_negative=(~mask_positive)
        parameter[mask_positive]=np.ceil((parameter[mask_positive]-eb)/(2*eb))
        parameter[mask_negative]=np.floor((parameter[mask_negative]+eb)/(2*eb))
        return parameter.astype(np.int32),2*eb*parameter
    rel_eb=1e-3
    abs_eb=rel_eb*(pressure.max()-pressure.min())
    parameter_eb=2*rel_eb*0.5
    f_qb,f_quantized=quantize_ndarray(f.copy(),parameter_eb)
    f=f_quantized

pressure2=np.einsum("ijk,kc->ijc",f,L)
pressure2=np.einsum("ijc,jb->ibc",pressure2,L)
pressure2=np.einsum("ibc,ia->abc",pressure2,L)

if lossy:
    hf_bs=b""
    pressure_qb,pressure_quantized=quantize_ndarray(pressure[0,0]-pressure2,abs_eb)
    pressure_quantized+=pressure2
    ######## Huffman f ########
    f_values,f_counts=np.unique(f_qb,return_counts=True)
    f_freq=dict(zip(f_values,f_counts))
    hf_bs+=struct.pack("<i",len(f_values))
    max_f=np.abs(f_values).max().item()
    if max_f<=2**15-1:
        hf_bs+=struct.pack("<b",15)
        hf_bs+=struct.pack(f"<{len(f_values)}h",*f_values.tolist())
    else:
        hf_bs+=struct.pack("<b",31)
        hf_bs+=struct.pack(f"<{len(f_values)}i",*f_values.tolist())
    hf_bs+=struct.pack(f"<{len(f_values)}i",*f_counts.tolist())
    f_codec=HuffmanCodec.from_frequencies(f_freq)
    temp_bs=f_codec.encode(f_qb.flatten().tolist())
    hf_bs+=struct.pack("<i",len(temp_bs))
    hf_bs+=temp_bs
    print(len(hf_bs))
    ######## Huffman qb ########
    qb_values,qb_counts=np.unique(pressure_qb,return_counts=True)
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
    temp_bs=qb_codec.encode(pressure_qb.flatten().tolist())
    hf_bs+=struct.pack("<i",len(temp_bs))
    hf_bs+=temp_bs
    print(len(hf_bs))
    ######## Zstd ########
    cctx=zstd.ZstdCompressor()
    zstd_bs=cctx.compress(hf_bs)
    print(len(zstd_bs))
    pressure2=pressure_quantized

if lossy:
    print(f"abs_eb= {abs_eb}")
print(f"max_err= {np.abs(pressure2-pressure[0,0]).max()}")
print(f"PSNR= {10*np.log10((pressure[0,0].max()-pressure[0,0].min())**2/np.mean((pressure2-pressure[0,0])**2))} dB")