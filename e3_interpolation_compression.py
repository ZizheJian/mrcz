import yaml,struct,os,time,math
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from read_nek5000_header import nek5000_header_c
from py_funcs.interpolation_compress_element import interpolation_compress_element
from py_funcs.predict_and_overwrite_linear import predict_and_overwrite_linear_dim0,predict_and_overwrite_linear_dim1,predict_and_overwrite_linear_dim2
from py_funcs.predict_and_overwrite_linear import predict_and_overwrite_linear_dim0_v2,predict_and_overwrite_linear_dim1_v2,predict_and_overwrite_linear_dim2_v2

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
# rel_eb=1e-5
abs_eb=1e-8

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

try:
    abs_eb
except:
    abs_eb=rel_eb*(header.data_range["P"][1]-header.data_range["P"][0]) # type: ignore
else:
    rel_eb=abs_eb/(header.data_range["P"][1]-header.data_range["P"][0])
max_err=0
mse=0
print(f"presure_min= {header.data_range['P'][0]}, presure_max= {header.data_range['P'][1]}")
print(f"rel_eb= {rel_eb}")
print(f"abs_eb= {abs_eb}")

x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header.element_shape[0]-1).deriv().roots(),[1]))

def compress_chunk(chunk):
    qb_list=[]
    max_err=0
    mse=0
    for element_id in range(chunk.shape[1]):
        element=chunk[0,element_id]
        qb=np.array([],dtype=np.int32)
        lvs=math.ceil(np.log2(element.shape[0]-1))
        for lv in range(lvs):
            if lv==0:
                previous_indices=[0,9]
            else:
                previous_indices=next_indices
            new_indices=[(previous_indices[i]+previous_indices[i+1])//2 if previous_indices[i+1]-previous_indices[i]>1 else None for i in range(len(previous_indices)-1)]
            next_indices=sorted(previous_indices+[i for i in new_indices if i is not None])
            local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim2(element,x,abs_eb,previous_indices,new_indices,next_indices)
            qb=np.concatenate((qb,local_qb))
            max_err=max(max_err,local_max_err)
            mse+=local_mse
            local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim1(element,x,abs_eb,previous_indices,new_indices,next_indices)
            qb=np.concatenate((qb,local_qb))
            max_err=max(max_err,local_max_err)
            mse+=local_mse
            local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim0(element,x,abs_eb,previous_indices,new_indices,next_indices)
            qb=np.concatenate((qb,local_qb))
            max_err=max(max_err,local_max_err)
            mse+=local_mse
        qb_list.append(qb)
    return qb_list,max_err,mse

def compress_chunk_v2(chunk):
    qb_list=[]
    max_err=0
    mse=0
    for element_id in range(chunk.shape[1]):
        element=chunk[0,element_id]
        qb=np.zeros_like(element,dtype=np.int32)
        lvs=math.ceil(np.log2(element.shape[0]-1))
        for lv in range(lvs):
            if lv==0:
                previous_indices=[0,9]
            else:
                previous_indices=next_indices
            new_indices=[(previous_indices[i]+previous_indices[i+1])//2 if previous_indices[i+1]-previous_indices[i]>1 else None for i in range(len(previous_indices)-1)]
            next_indices=sorted(previous_indices+[i for i in new_indices if i is not None])
            local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim2_v2(element,x,abs_eb,previous_indices,new_indices,next_indices)
            if np.any(np.logical_and(qb!=0,local_qb!=0)):
                raise ValueError("Overlap in qb!")
            qb+=local_qb
            max_err=max(max_err,local_max_err)
            mse+=local_mse
            local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim1_v2(element,x,abs_eb,previous_indices,new_indices,next_indices)
            if np.any(np.logical_and(qb!=0,local_qb!=0)):
                raise ValueError("Overlap in qb!")
            qb+=local_qb
            max_err=max(max_err,local_max_err)
            mse+=local_mse
            local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim0_v2(element,x,abs_eb,previous_indices,new_indices,next_indices)
            if np.any(np.logical_and(qb!=0,local_qb!=0)):
                raise ValueError("Overlap in qb!")
            qb+=local_qb
            max_err=max(max_err,local_max_err)
            mse+=local_mse
        qb_list.append(qb)
    return qb_list,max_err,mse

# def compress_chunk_sz3(chunk):
#     qb_list=[]
#     max_err=0
#     mse=0
#     for element_id in range(chunk.shape[1]):
#         element=chunk[0,element_id]
#         qb=np.array([],dtype=np.int32)
#         lvs=math.ceil(np.log2(element.shape[0]-1))
#         for lv in range(lvs):
#             if lv==0:
#                 previous_indices=[0,9]
#             else:
#                 previous_indices=next_indices
#             new_indices=[(previous_indices[i]+previous_indices[i+1])//2 if previous_indices[i+1]-previous_indices[i]>1 else None for i in range(len(previous_indices)-1)]
#             next_indices=sorted(previous_indices+[i for i in new_indices if i is not None])
#             local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim2(element,x,abs_eb,previous_indices,new_indices,next_indices)
#             qb=np.concatenate((qb,local_qb))
#             max_err=max(max_err,local_max_err)
#             mse+=local_mse
#             local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim1(element,x,abs_eb,previous_indices,new_indices,next_indices)
#             qb=np.concatenate((qb,local_qb))
#             max_err=max(max_err,local_max_err)
#             mse+=local_mse
#             local_qb,local_max_err,local_mse=predict_and_overwrite_linear_dim0(element,x,abs_eb,previous_indices,new_indices,next_indices)
#             qb=np.concatenate((qb,local_qb))
#             max_err=max(max_err,local_max_err)
#             mse+=local_mse
#         qb_list.append(qb)
#     return qb_list,max_err,mse

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
    total_lossy_time=0
    for block in blocks:
        print(f"block {block['block_id']+1}/{block_num}")
        block_element_num=os.path.getsize(block["pressure_path"])//(header.word_size*header.element_size)
        with open(block["pressure_path"],"rb") as f_pressure:
            cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=block_element_num*header.element_size).reshape(1,-1,*header.element_shape)
        start_element_id=[(block_element_num*i)//thread_num for i in range(thread_num)]
        end_element_id=[(block_element_num*(i+1))//thread_num for i in range(thread_num)]
        exe_input=[cache[:,start_element_id[i]:end_element_id[i]] for i in range(thread_num)]
        start_time=time.time()
        results=list(exe.map(compress_chunk_v2,exe_input))
        total_lossy_time+=time.time()-start_time
        qb_list.append(([],))
        for local_qb_list,local_max_err,local_mse in results:
            qb_list[-1][0].extend(local_qb_list)
            if local_max_err>max_err:
                max_err=local_max_err
            mse+=local_mse
    print(f"有损时间: {total_lossy_time} s")
    
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