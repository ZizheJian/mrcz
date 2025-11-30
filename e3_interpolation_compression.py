import yaml,struct,os,time
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from read_nek5000_header import nek5000_header_c
from py_funcs.interpolation_compress_element import interpolation_compress_element

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

# rel_eb=1e-5
abs_eb=1e-5
try:
    rel_eb # type: ignore
except:
    rel_eb=abs_eb/(header.data_range["P"][1]-header.data_range["P"][0])
else:
    abs_eb=rel_eb*(header.data_range["P"][1]-header.data_range["P"][0]) # type: ignore
max_err=0
mse=0
print(f"presure_min={header.data_range['P'][0]}, presure_max={header.data_range['P'][1]}")
print(f"rel_eb={rel_eb}")
print(f"abs_eb={abs_eb}")

x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header.element_shape[0]-1).deriv().roots(),[1]))

def compress_chunk(args):
    block,start_element_id,end_element_id=args
    with open(block["pressure_path"],"rb") as f_pressure:
        f_pressure.seek(start_element_id*header.element_size*header.word_size)
        cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=(end_element_id-start_element_id)*header.element_size)
        cache=cache.reshape((1,-1)+header.element_shape)
    qb_list=[]
    max_err=0
    mse=0
    for element_id in range(end_element_id-start_element_id):
        element=cache[0,element_id]
        qb,local_max_err,local_mse=interpolation_compress_element(element,header.element_shape,x,abs_eb)
        qb_list.append(qb)
        max_err=max(max_err,local_max_err)
        mse+=local_mse
    return qb_list,max_err,mse

def lossless_encode_block(qb_list):
    qb=np.array(qb_list).flatten()
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
    ######## Zstd ########
    cctx=zstd.ZstdCompressor()
    zstd_bs=cctx.compress(hf_bs)
    return len(zstd_bs)

qb_list=[]
with ProcessPoolExecutor(max_workers=thread_num) as exe:
    start_time=time.time()
    for block in blocks:
        print(f"block {block['block_id']+1}/{block_num}")
        block_element_num=os.path.getsize(block["pressure_path"])//(header.word_size*header.element_size)
        exe_input=[(block,(block_element_num*i)//thread_num,(block_element_num*(i+1))//thread_num)
                     for i in range(thread_num)]
        results=list(exe.map(compress_chunk,exe_input))
        qb_list.append(([],))
        for local_qb_list,local_max_err,local_mse in results:
            qb_list[-1][0].extend(local_qb_list)
            if local_max_err>max_err:
                max_err=local_max_err
            mse+=local_mse
    print(f"有损时间: {time.time()-start_time} s")
    
    start_time=time.time()
    results=list(exe.map(lossless_encode_block,qb_list))
    zstd_bs_len=sum(results)
    print(f"无损时间: {time.time()-start_time} s")

print(f"Zstd编码长度: {zstd_bs_len}")
print(f"原始数据大小: {header.file_element_num*header.element_size*header.word_size} bytes")
print(f"压缩比: {header.file_element_num*header.element_size*header.word_size/zstd_bs_len}")

print(f"max_err= {max_err}({max_err*100/abs_eb}%)")
mse/=(header.file_element_num*header.element_size)
psnr=10*np.log10((header.data_range["P"][1]-header.data_range["P"][0])**2/mse)
print(f"MSE= {mse} , PSNR= {psnr} dB")