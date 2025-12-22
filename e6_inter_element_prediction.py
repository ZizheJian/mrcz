import yaml,struct,os,time,math
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from read_nek5000_header import nek5000_header_c
from py_funcs.quantize import quantize_ndarray1,quantize_ndarray2,quantize_float1,quantize_float2

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
# rel_eb=1e-5
abs_eb=1e-6

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

def compress_line(line):
    qb=np.zeros_like(line,dtype=np.int32)
    new_line=line.copy()
    max_err=0
    mse=0
    lvs=math.ceil(np.log2(line.shape[0]))
    qb[0],new_line[0]=quantize_float1(line[0],abs_eb)
    qb_index=1
    for lv in range(lvs):
        step_size=2**(lvs-lv-1)
        end_index=(line.shape[0]+step_size-1)//step_size
        if end_index==1:
            continue
        elif end_index==2:
            qb[qb_index],new_line[step_size]=quantize_float2(line[0],line[step_size],abs_eb)
            qb_index+=1
        elif end_index>=3:
            pred_num=(end_index-1)//2
            qb[qb_index:qb_index+pred_num],new_line[step_size:step_size*(2*pred_num+1):2*step_size]= \
                quantize_ndarray2((line[0:step_size*2*pred_num:2*step_size]+line[step_size*2:step_size*(2*pred_num+2):2*step_size])/2,
                                  line[step_size:step_size*(2*pred_num+1):2*step_size],abs_eb)
            qb_index+=pred_num
            # for i in range(1,end_index-1,2):
            #     qb[qb_index],new_line[i*step_size]=quantize_float2((line[(i-1)*step_size]+line[(i+1)*step_size])/2,line[i*step_size],abs_eb)
            #     qb_index+=1
            if end_index%2==0:
                qb[qb_index],new_line[(end_index-1)*step_size]=quantize_float2(line[(end_index-2)*step_size],line[(end_index-1)*step_size],abs_eb)
                qb_index+=1
    max_err=np.abs(new_line-line).max()
    mse=np.sum((new_line-line)**2)
    return qb,max_err,mse

def lossless_encode_block(qb_list):
    qb=np.array(qb_list).flatten()
    print(qb.shape)
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
        print(f"Block {block['block_id']+1}/{block_num}")
        block_element_num=os.path.getsize(block["pressure_path"])//(header.word_size*header.element_size)
        with open(block["pressure_path"],"rb") as f_pressure:
            cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=block_element_num*header.element_size).reshape(1,-1,*header.element_shape).astype(np.float64)
        cache=cache.reshape(1,cache.shape[1],-1)
        exe_input=[cache[0,:,i] for i in range(cache.shape[2])]
        start_time=time.time()
        results=list(exe.map(compress_line,exe_input))
        total_lossy_time+=time.time()-start_time
        qb_list.append([])
        for local_qb_ndarray,local_max_err,local_mse in results:
            qb_list[-1].append(local_qb_ndarray)
            if local_max_err>max_err:
                max_err=local_max_err
            mse+=local_mse
    print(f"有损时间: {total_lossy_time} s")

    # for block_id,(f_qb_list,qb_list) in enumerate(qb_list):
    #     print(f"Block {block_id}")
    #     f_qb=np.array(f_qb_list)
    #     f_qb_0_freq=np.sum(f_qb==0,axis=0)/f_qb.shape[0]

    #     fig,axes=plt.subplots(2,5,figsize=(50,20))
    #     axes=axes.flatten()
    #     ims=[]
    #     for i in range(10):
    #         im=axes[i].imshow(f_qb_0_freq[i,:,:],origin="lower",cmap="gray",vmin=0,vmax=1)
    #         ims.append(im)
    #         axes[i].axis("off")
    #     cbar=fig.colorbar(ims[0],ax=axes,fraction=0.05,shrink=0.9,aspect=20,pad=0.05)
    #     cbar.set_ticks([0,1])
    #     cbar.set_ticklabels(["No Zero","All Zero"],fontsize=60)
    #     plt.savefig(f"fig/0_freq_{block_id}.png")
    #     plt.close()

    #     x=y=z=np.arange(header.element_shape[0])
    #     X,Y,Z=np.meshgrid(x,y,z,indexing="ij")
    #     fig=go.Figure()
    #     fig.add_trace(go.Volume(x=X.flatten(),y=Y.flatten(),z=Z.flatten(),value=f_qb_0_freq.flatten(),isomin=0,isomax=1,
    #                             opacity=0.1,surface_count=40,colorscale="Greys",showscale=True,
    #                             colorbar=dict(tickmode="array",tickvals=[0,1],ticktext=["No Zero","All Zero"])))
    #     fig.update_layout(scene=dict(aspectmode="cube"),width=800,height=800)
    #     fig.write_html(f"fig/0_freq_{block_id}.html")
    
    start_time=time.time()
    results=list(exe.map(lossless_encode_block,qb_list))
    print(f"无损时间: {time.time()-start_time} s")
    zstd_len=sum(results)

print(f"Zstd编码长度: {zstd_len}")
print(f"原始数据大小: {header.file_element_num*header.element_size*header.word_size} bytes")
print(f"压缩比: {header.file_element_num*header.element_size*header.word_size/zstd_len}")

print(f"max_err= {max_err}({max_err*100/abs_eb}%)")
mse/=(header.file_element_num*header.element_size)
psnr=10*np.log10((header.data_range["P"][1]-header.data_range["P"][0])**2/mse)
print(f"MSE= {mse} , PSNR= {psnr} dB")