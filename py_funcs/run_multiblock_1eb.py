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
from py_funcs.global_param_class import global_param_class

def compress_chunk(args):
    chunk,param_eb,abs_eb,zero_freq_bound,L,neg_L,element_base_weight=args
    f_qb_list=[]
    qb_list=[]
    max_err=0
    mse=0
    for element_id in range(chunk.shape[1]):
        element=chunk[0,element_id]
        element_base=element[0,0,0]*element_base_weight[0]+element[0,0,-1]*element_base_weight[1]+\
                    element[0,-1,0]*element_base_weight[2]+element[0,-1,-1]*element_base_weight[3]+\
                    element[-1,0,0]*element_base_weight[4]+element[-1,0,-1]*element_base_weight[5]+\
                    element[-1,-1,0]*element_base_weight[6]+element[-1,-1,-1]*element_base_weight[7]
        element-=element_base
        f=np.einsum("abc,ai->ibc",element,neg_L)
        f=np.einsum("ibc,bj->ijc",f,neg_L)
        f=np.einsum("ijc,ck->ijk",f,neg_L)
        f_qb,f_quantized=quantize_ndarray3(f,param_eb,zero_freq_bound)
        f_qb_list.append(f_qb)
        new_element=np.einsum("ijk,kc->ijc",f_quantized,L)
        new_element=np.einsum("ijc,jb->ibc",new_element,L)
        new_element=np.einsum("ibc,ia->abc",new_element,L)
        qb,new_element=quantize_ndarray2(new_element,element,abs_eb)
        qb_list.append(qb)
        max_err=max(max_err,np.abs(new_element-element).max())
        mse+=np.sum((new_element-element)**2)
    return f_qb_list,qb_list,max_err,mse

def lossless_encode_block(args):
    head_bs=b""
    f_qb_postprocessing_method=args[0]
    if f_qb_postprocessing_method=="None":
        f_qb_list,qb_list=args[1:]
        f_qb=np.array(f_qb_list).transpose(1,2,3,0).flatten()
        qb=np.array(qb_list).flatten()
    else:
        f_qb_list,qb_list,bitmap=args[1:]
        if np.all(bitmap):
            head_bs+=struct.pack("B",0)
            return 1
        elif np.all(bitmap==False):
            head_bs+=struct.pack("B",1)
            f_qb=np.array(f_qb_list).transpose(1,2,3,0).flatten()
        else:
            head_bs+=struct.pack("B",2)
            f_qb=np.array(f_qb_list)
            f_qb=f_qb[:,bitmap==False].transpose(1,0).flatten()
        qb=np.array(qb_list).flatten()
    ######## Huffman f ########
    f_hf_bs=b""
    f_values,f_counts=np.unique(f_qb,return_counts=True)
    f_freq=dict(zip(f_values,f_counts))
    f_hf_bs+=struct.pack("<Q",len(f_values))
    max_f=np.abs(f_values).max().item()
    if max_f<=2**15-1:
        f_hf_bs+=struct.pack("<b",15)
        f_hf_bs+=struct.pack(f"<{len(f_values)}h",*f_values.tolist())
    elif max_f<=2**31-1:
        f_hf_bs+=struct.pack("<b",31)
        f_hf_bs+=struct.pack(f"<{len(f_values)}i",*f_values.tolist())
    else:
        f_hf_bs+=struct.pack("<b",63)
        f_hf_bs+=struct.pack(f"<{len(f_values)}q",*f_values.tolist())
    f_hf_bs+=struct.pack(f"<{len(f_values)}Q",*f_counts.tolist())
    f_codec=HuffmanCodec.from_frequencies(f_freq)
    temp_bs=f_codec.encode(f_qb.flatten().tolist())
    f_hf_bs+=struct.pack("<Q",len(temp_bs))
    f_hf_bs+=temp_bs
    ######## Huffman qb ########
    hf_bs=b""
    qb_values,qb_counts=np.unique(qb,return_counts=True)
    qb_freq=dict(zip(qb_values,qb_counts))
    hf_bs+=struct.pack("<Q",len(qb_values))
    max_qb=np.abs(qb_values).max().item()
    if max_qb<=2**15-1:
        hf_bs+=struct.pack("<b",15)
        hf_bs+=struct.pack(f"<{len(qb_values)}h",*qb_values.tolist())
    elif max_qb<=2**31-1:
        hf_bs+=struct.pack("<b",31)
        hf_bs+=struct.pack(f"<{len(qb_values)}i",*qb_values.tolist())
    else:
        hf_bs+=struct.pack("<b",63)
        hf_bs+=struct.pack(f"<{len(qb_values)}q",*qb_values.tolist())
    hf_bs+=struct.pack(f"<{len(qb_values)}Q",*qb_counts.tolist())
    qb_codec=HuffmanCodec.from_frequencies(qb_freq)
    temp_bs=qb_codec.encode(qb.flatten().tolist())
    hf_bs+=struct.pack("<Q",len(temp_bs))
    hf_bs+=temp_bs
    ######## Zstd ########
    cctx=zstd.ZstdCompressor()
    zstd_bs=cctx.compress(head_bs+f_hf_bs+hf_bs)
    return len(zstd_bs)

def run_multiblock_1eb(global_param:global_param_class,block_list:list[int],abs_eb:float)->tuple[list[float],list[float],list[float],list[float]]:
    param_eb_scale=global_param.param["param_eb_scale"]
    zero_freq_bound=global_param.param["zero_freq_bound"]
    L=global_param.param["L"]
    neg_L=global_param.param["neg_L"]
    element_base_weight=global_param.param["element_base_weight"]
    header=global_param.header
    thread_num=global_param.thread_num
    f_qb_postprocessing_method=global_param.param["f_qb_postprocessing_method"]

    param_eb=abs_eb*param_eb_scale
    max_err=0

    qb_list=[]
    all_zero_bitmap_list=[] # True = 在element内某个位置上，沿着element_id轴的f_qb全为0
    mse_list=[]
    with ProcessPoolExecutor(max_workers=global_param.thread_num) as exe:
        for block_id in block_list:
            block=global_param.blocks[block_id]
            block_element_num=block["element_num"]
            with open(block["pressure_path"],"rb") as f_pressure:
                cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=block_element_num*header.element_size).reshape(1,-1,*header.element_shape).astype(np.float64)
            start_element_id=[(block_element_num*i)//thread_num for i in range(thread_num)]
            end_element_id=[(block_element_num*(i+1))//thread_num for i in range(thread_num)]
            other_input_param=(param_eb,abs_eb,zero_freq_bound,L,neg_L,element_base_weight)
            exe_input=[(cache[:,start_element_id[i]:end_element_id[i]],)+other_input_param for i in range(thread_num)]
            results=list(exe.map(compress_chunk,exe_input))
            qb_list.append(([],[]))
            mse_list.append(0)
            for local_f_qb_list,local_qb_list,local_max_err,local_mse in results:
                qb_list[-1][0].extend(local_f_qb_list)
                qb_list[-1][1].extend(local_qb_list)
                if local_max_err>max_err:
                    max_err=local_max_err
                mse_list[-1]+=local_mse

        if f_qb_postprocessing_method=="bitmap":
            for local_f_qb_list,local_qb_list in qb_list:
                f_qb=np.array(local_f_qb_list)
                all_zero_bitmap_list.append(np.all(f_qb==0,axis=0))
        if f_qb_postprocessing_method=="None":
            exe_input=[(f_qb_postprocessing_method,qb_list[block_id][0],qb_list[block_id][1]) for block_id in range(len(block_list))]
            results=list(exe.map(lossless_encode_block,exe_input))
        else:
            exe_input=[(f_qb_postprocessing_method,qb_list[block_id][0],qb_list[block_id][1],all_zero_bitmap_list[0]) for block_id in range(len(block_list))]
            results=list(exe.map(lossless_encode_block,exe_input))
        zstd_len_list=list(results)

    if max_err>abs_eb*1.001:
        raise ValueError(f"Compression error {max_err} exceeds the absolute error bound {abs_eb}!")
    br_list=[zstd_len_list[i]/(global_param.blocks[block_list[i]]["element_num"]*header.element_size) for i in range(len(block_list))]
    cr_list=[header.word_size/br for br in br_list]
    mse_list=[mse_list[i]/(global_param.blocks[block_list[i]]["element_num"]*header.element_size) for i in range(len(block_list))]
    psnr_list=[10*np.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse) for mse in mse_list]
    return br_list,cr_list,mse_list,psnr_list