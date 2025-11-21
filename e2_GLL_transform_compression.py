import yaml,struct,os,random
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from read_nek5000_header import nek5000_header_c

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

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

# rel_eb=1e-5
# abs_eb=(header.data_range["P"][1]-header.data_range["P"][0])*rel_eb
abs_eb=1e-5
rel_eb=abs_eb/(header.data_range["P"][1]-header.data_range["P"][0])
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

def compress_chunk(args):
    block,start_element_id,end_element_id=args
    with open(block["pressure_path"],"rb") as f_pressure:
        f_pressure.seek(start_element_id*header.element_size*header.word_size)
        cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=(end_element_id-start_element_id)*header.element_size)
        cache=cache.reshape((1,-1)+header.element_shape)
    f_qb_list=[]
    qb_list=[]
    max_err=0
    mse=0
    for element_id in range(end_element_id-start_element_id):
        element=cache[0,element_id]
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
    return f_qb_list,qb_list,max_err,mse

f_qb_list=[]
qb_list=[]
selected_element_num=0
with ProcessPoolExecutor(max_workers=thread_num) as exe:
    for block in blocks:
        # if random.random()>0.1:
        #     continue
        print(f"block {block['block_id']+1}/{block_num}")
        block_element_num=os.path.getsize(block["pressure_path"])//(header.word_size*header.element_size)
        selected_element_num+=block_element_num
        exe_input=[(block,(block_element_num*i)//thread_num,(block_element_num*(i+1))//thread_num)
                     for i in range(thread_num)]
        results=list(exe.map(compress_chunk,exe_input))
        for local_f_qb_list,local_qb_list,local_max_err,local_mse in results:
            f_qb_list.extend(local_f_qb_list)
            qb_list.extend(local_qb_list)
            if local_max_err>max_err:
                max_err=local_max_err
            mse+=local_mse
f_qb=np.array(f_qb_list).flatten()
f_qb=f_qb[f_qb!=0]
qb=np.array(qb_list).flatten()
print(f_qb.shape,qb.shape)

######## Huffman f ########
f_hf_bs=b""
f_values,f_counts=np.unique(f_qb,return_counts=True)
f_freq=dict(zip(f_values,f_counts))
f_hf_bs+=struct.pack("<Q",len(f_values))
max_f=np.abs(f_values).max().item()
if max_f<=2**15-1:
    f_hf_bs+=struct.pack("<b",15)
    f_hf_bs+=struct.pack(f"<{len(f_values)}h",*f_values.tolist())
else:
    f_hf_bs+=struct.pack("<b",31)
    f_hf_bs+=struct.pack(f"<{len(f_values)}i",*f_values.tolist())
f_hf_bs+=struct.pack(f"<{len(f_values)}Q",*f_counts.tolist())
f_codec=HuffmanCodec.from_frequencies(f_freq)
temp_bs=f_codec.encode(f_qb.flatten().tolist())
f_hf_bs+=struct.pack("<Q",len(temp_bs))
f_hf_bs+=temp_bs
print(f"频率编码长度: {len(f_hf_bs)}")
######## Huffman qb ########
hf_bs=b""
qb_values,qb_counts=np.unique(qb,return_counts=True)
qb_freq=dict(zip(qb_values,qb_counts))
hf_bs+=struct.pack("<Q",len(qb_values))
max_qb=np.abs(qb_values).max().item()
if max_qb<=2**15-1:
    hf_bs+=struct.pack("<b",15)
    hf_bs+=struct.pack(f"<{len(qb_values)}h",*qb_values.tolist())
else:
    hf_bs+=struct.pack("<b",31)
    hf_bs+=struct.pack(f"<{len(qb_values)}i",*qb_values.tolist())
hf_bs+=struct.pack(f"<{len(qb_values)}Q",*qb_counts.tolist())
qb_codec=HuffmanCodec.from_frequencies(qb_freq)
temp_bs=qb_codec.encode(qb.flatten().tolist())
hf_bs+=struct.pack("<Q",len(temp_bs))
hf_bs+=temp_bs
print(f"qb编码长度: {len(hf_bs)}")
######## Zstd ########
cctx=zstd.ZstdCompressor()
zstd_bs=cctx.compress(f_hf_bs+hf_bs)
print(f"Zstd编码长度: {len(zstd_bs)}")
print(f"原始数据大小: {selected_element_num*header.element_size*header.word_size} bytes")
print(f"压缩比: {selected_element_num*header.element_size*header.word_size/len(zstd_bs)}")
# print(f"原始数据大小: {header.file_element_num*header.element_size*header.word_size} bytes")
# print(f"压缩比: {header.file_element_num*header.element_size*header.word_size/len(zstd_bs)}")

print(f"max_err= {max_err}({max_err*100/abs_eb}%)")
mse/=(header.file_element_num*header.element_size)
psnr=10*np.log10((header.data_range["P"][1]-header.data_range["P"][0])**2/mse)
print(f"MSE= {mse} , PSNR= {psnr} dB")